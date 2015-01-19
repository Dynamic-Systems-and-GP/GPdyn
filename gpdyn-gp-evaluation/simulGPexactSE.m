function [m, s2, mu, sig2] = simulGPexactSE(hyp, inf, mean, cov, lik, input, target, test, lag)
% Simulation of the GP model, where the output variance is propagated using
% analytical approximation
%
%% Syntax
%  [m, s2, mu, sig2] = simulGPexactSE(hyp, inf, mean, cov, lik, input, target, test, lag)
%
%% Description
% See A. Girard, Approximate Methods for Propagation of
% Uncertainty with Gaussian Process Models, PhD thesis, 2004. 
% Simulation of the GP model, where the output variance is propagated using
% analytical approximation. It can be used only with Gaussian 
% covariance function and Gaussian likelihood function with white noise
% model (covSEard and likGauss).  
% Uses routine gpExactSEard. 
% 
% Inputs: 
% *  hyp      ... struct of optimized hyperparameters 
% *  inf      ... function specifying the inference method 
% *  mean     ... prior mean function
% *  cov      ... specified covariance function, see help covFun for more info 
% *  lik      ... likelihood function
% *  input    ... input part of the training data,  NxD matrix
% *  target   ... output part of the training data (ie. target), Nx1 vector 
% *  test     ... input matrix for simulation, kxD vector, see
%                 construct.m for more info 
% *  lag      ... the order of the model (number of used lagged outputs) 
% 
% Outputs: 
% *  m     ... predictive mean when propagating the uncertainty 
% *  s2    ... predictive variance when propagating the uncertainty 
% *  mu    ... predictive mean using "naive" approach (doesn't propagate
%              the uncertainty)
% *  sig2  ... predictive variance using "naive" approach
% 
% See also:
% gpExactSEard, simulGPtaylorSE, simulGPnaive
%
% Examples:
% demo_example_gp_simulation.m
%

%% 
% Written by K. Azman, 2007, based on the software of J. Quinonero-Candela 
% and A. Girard. 

fun_name = 'simulGPexactSE'; 


[n, D] = size(input);
[nn, D] = size(test);

% input validation
[ is_valid, hyp, inf, mean, cov, lik, msg ] = validate( hyp, inf, mean, cov, lik, D);

if ~isequal(cov,{@covSEard}) 
    error(strcat([fun_name,': function can only be called with the', ...
        ' covariance function ''covSEard'' '])); 
end

if ~isequal(lik,{@likGauss}) 
    error(strcat([fun_name,': function can only be called with the', ...
        ' likelihood function ''likGauss'', where hyp.lik parameter is log(sn)'])); 
end 



% if ~(isequal(cov{1},'covSum') &  isequal(cov{2}{1},'covSEard') & ...
%         isequal(cov{2}{2},'covNoise'))
%     error(strcat([fun_name,': function can be called only with the sum', ...
%         ' of covariance functions ''covSEard'' and ''covNoise'' '])); 
% end 

X=[-2*hyp.cov(1:end-1);2*hyp.cov(end);2*hyp.lik]; % adapt hyperparameters to local format
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
muXp = test(1,:); % propagate uncertainty
SigX = zeros(D,D); 
muX = test(1,:);  % naive approach

if (nargin==9) % in the case of control inputs
    [mtmp, s2tmp] = gpExactSEard(hyp, inf, mean,cov,lik, invQ, input, target, muXp, SigX, lag);
else
    [mtmp, s2tmp] = gpExactSEard(hyp, inf,mean,cov,lik, invQ, input, target, muXp, SigX);
end
m(1) = mtmp; 

s2(1) = s2tmp; 

[mutmp, sig2tmp] = gpExactSEard(hyp, inf, mean,cov,lik, invQ, input, target, muX);
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
        if (nargin==9) % in the case of control inputs
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
    
    if (nargin==9) % control inputs to take into account
        SigX(1:lag-1,1:lag-1) = SigX(2:lag,2:lag);
        SigX(lag,lag) = s2(k-1);
        SigX(1:lag-1,lag) = covXY(2:lag);
        SigX(lag,1:lag-1) = covXY(2:lag)';

        
        muXp = [muXp(2:lag) m(k-1) test(k,lag+1:end)];    
        muX = [muX(2:lag) mu(k-1) test(k,lag+1:end)];
        [mtmp, s2tmp] = gpExactSEard(hyp, inf, mean,cov,lik, invQ, input, target, muXp, SigX, lag);
    else
        SigX(1:D-1,1:D-1) = SigX(2:D,2:D);
        SigX(D,D) = s2(k-1);
        SigX(1:D-1,D) = covXY(2:D); 
        SigX(D,1:D-1) = covXY(2:D)'; 

        muXp = [muXp(2:D) m(k-1)];    
        muX = [muX(2:D) mu(k-1)];
       [mtmp, s2tmp] = gpExactSEard(hyp, inf, mean,cov,lik, invQ, input, target, muXp, SigX);
    end
         
    m(k) = mtmp; 
    
    if s2tmp < 0
        
        warning(strcat([fun_name,': negative variance, step:',num2str(k)]))
   
    end

    s2(k) = s2tmp;
    
    
    [mutmp, sig2tmp] = gpExactSEard(hyp, inf,mean, cov,lik, invQ, input, target, muX);

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
% output noise variance added in the end
 s2 = s2 + vy;
% output noise variance added in the end 
 sig2 = sig2 + vy; 



return; 

