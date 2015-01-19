%% simul00exact
function [mu, sig2, m, s2] = simul00exact(logtheta, covfunc, input, target, xt, lag)

%% Syntax
%  [mu, sig2, m, s2] = simul00exact(logtheta, covfunc, input, target, xt, lag)

%% Description
% Simulation of the GP model, where the output variance is propagated using
% analytical approximation, see A. Girard, Approximate Methods for Propagation of
% Uncertainty with Gaussian Process Models, PhD thesis, 2004. 
% Notes:
% Currently it can be used only with Gaussian covariance function and
% with white noise model (sum of covSEard and covNoise). 

% Uses routine gpr_SEard_exact. 
% 
% Inputs: 
% loghteta .. optimized hyperparameters 
% covfunc .. specified covariance function, see help covFun for more info 
% input .. input part of the training data,  NxD matrix
% target .. output part of the training data (ie. target), Nx1 vector 
% xt .. input matrix for simulation, kxD vector, see
%   construct_simul_input.m for more info 
% lag .. the order of the model (number of used lagged outputs) 
% Outputs: 
% mu        predictive mean using "naive" approach (doesn't propagate the
%           uncertainty)
% sig2      predictive variance using "naive" approach
% m         predictive mean when propagating the uncertainty 
% s2        predictive variance when propagating the uncertainty 
% The form of the covariance function is
%    C(x^p,x^q) = v1*exp( -0.5 * sum_{d=1..D} w_d * (x^p_d - x^q_d)^2 ) + v1
%    \delta_{pq}
% (computed using cov2.m)
%
% Based on the work of J. Quinonero-Candela and A. Girard. 

%% Examples
% demo_example_gp_simulation.m

%% See Also
% SIMUL00NAIVE, GPR_SEARD_EXACT

fun_name = 'simul00exact'; 


if ~(isequal(covfunc{1},'covSum') &  isequal(covfunc{2}{1},'covSEard') & ...
        isequal(covfunc{2}{2},'covNoise'))
    error(strcat([fun_name,': function can be called only with the sum', ...
        ' of covariance functions ''covSEard'' and ''covNoise'' '])); 
end 




[n, D] = size(input);
[nn, D] = size(xt);
X=[-2*logtheta(1:end-2);2*logtheta(end-1);2*logtheta(end)];
expX = exp(X);

vy = expX(end);

% training covariance matrix
Q = zeros(n,n);
for d = 1:D
  Q = Q + (repmat(input(:,d),1,n)-repmat(input(:,d)',n,1)).^2*expX(d);
end
Q = expX(D+1)*exp(-0.5*Q);

Q = Q + vy*eye(n);  % data cov. matrix: add noise

invQ = inv(Q);

beta = invQ*target;
b = expX(D+1); % v_1 really...

% 1st point
muXp = xt(1,:); % propagate uncertainty
SigX = zeros(D,D); 
muX = xt(1,:);  % naive approach

if (nargin==6) % in the case of control inputs
    [mtmp, s2tmp] = gpr_SEard_exact(logtheta, covfunc, invQ, input, target, muXp, SigX, lag);
else
    [mtmp, s2tmp] = gpr_SEard_exact(logtheta, covfunc, invQ, input, target, muXp, SigX);
end
m(1) = mtmp; 

s2(1) = s2tmp; 

[mutmp, sig2tmp] = gpr_SEard_exact(logtheta, covfunc, invQ, input, target, muX);
mu(1) = mutmp; 

sig2(1) = sig2tmp; 

for k=2:nn
    
    if (mod(k,50) == 0)
        disp(strcat([fun_name, ', step: ', int2str(k), '/', int2str(nn)]));
    end    
    
    % For the NEXT prediction...
    
    % cross-cov terms
    % little things to do before we start
    L = sum(diag(SigX)~=0); % number of stochastic dimensions (non zero vars)
    covXY = zeros(D,1);
    if L>0
        if (nargin==6) % in the case of control inputs
            rangeL = lag-L+1:lag;
         
        else
            rangeL = D-L+1:D;
       
        end
        rangeC = [1:rangeL(1)-1 rangeL(end)+1:D];
        
        SigXL = SigX(rangeL,rangeL);
        muXL = muXp(:,rangeL);
        muXC = muXp(:,rangeC);
        inputL = input(:,rangeL);
        inputC = input(:,rangeC);
        
        invLL = diag(expX(rangeL));
        invLC = diag(expX(rangeC));
        invSL = inv(SigXL);
        invC = (invLL+invSL);    
        invSmuX = invSL*muXL';
        t1 = muXL*invSmuX;
        c = invC\(invLL*inputL'+repmat(invSmuX,1,n));
        t2 = sum(inputL.*(inputL*invLL),2);
        t3 = sum(c.*(invC*c),1)';
        I = (1/sqrt(det(invLL*SigXL+eye(L))))*exp(-0.5*(t1+t2-t3));    
        
        CC = exp(-.5*(sum((inputC-repmat(muXC,n,1)).*((inputC-repmat(muXC,n,1))*invLC),2)));
        
        aux = m(k-1)*muXp';
        covXY(rangeC) = repmat(muXC',1,n)*(beta.*b.*(CC.*I));
        covXY(rangeL) = b*c*(beta.*(CC.*I));
        covXY = covXY - aux;
        covXY(rangeC) = zeros(length(rangeC),1);
  
        
    end

    % input covariance matrix and mean
    
    if (nargin==6) % control inputs to take into account
        SigX(1:lag-1,1:lag-1) = SigX(2:lag,2:lag);
        SigX(lag,lag) = s2(k-1);
        SigX(1:lag-1,lag) = covXY(2:lag);
        SigX(lag,1:lag-1) = covXY(2:lag)';

        
        muXp = [muXp(2:lag) m(k-1) xt(k,lag+1:end)];    
        muX = [muX(2:lag) mu(k-1) xt(k,lag+1:end)];
        [mtmp, s2tmp] = gpr_SEard_exact(logtheta, covfunc, invQ, input, target, muXp, SigX, lag);
    else
        SigX(1:D-1,1:D-1) = SigX(2:D,2:D);
        SigX(D,D) = s2(k-1);
        SigX(1:D-1,D) = covXY(2:D); 
        SigX(D,1:D-1) = covXY(2:D)'; 

        muXp = [muXp(2:D) m(k-1)];    
        muX = [muX(2:D) mu(k-1)];
       [mtmp, s2tmp] = gpr_SEard_exact(logtheta, covfunc, invQ, input, target, muXp, SigX);
    end
         
    m(k) = mtmp; 
    
    if s2tmp < 0
        
        warning(strcat([fun_name,': negative variance, step:',num2str(k)]))
   
    end

    s2(k) = s2tmp;
    
    
    [mutmp, sig2tmp] = gpr_SEard_exact(logtheta, covfunc, invQ, input, target, muX);

    mu(k) = mutmp;   
    
    if sig2tmp < 0
        
        warning(strcat([fun_name,': no propagation - negative variance, step:',...
            num2str(k)])); 
     
    end
    
    sig2(k) = sig2tmp; 
end   

% transform into coloumn vectors 
m = m'; 
s2 = s2'; 
mu = mu';
sig2 = sig2'; 

% in the case of negative variance, but not to cover it!
for i=1:length(s2)
    if s2(i)<0
        s2(i)=0;
    end
    if sig2(i)<0
        sig2(i)=0;
    end
end
% Azman, 16.10., output noise variance added in the end!!!! 
 s2 = s2 + vy; 
% Azman, 16.10., output noise variance added in the end!!!! 
 sig2 = sig2 + vy; 



return; 

