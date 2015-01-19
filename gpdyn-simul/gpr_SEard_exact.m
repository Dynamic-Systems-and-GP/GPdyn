%% gpr_SEard_exact
function [m, S2] = gpr_SEard_exact(logtheta, covfunc, invQ, input, target, muX, SigX, lag)

%% Syntax
%  [m, S2] = gpr_SEard_exact(logtheta, covfunc, invQ, input, target, muX, SigX, lag);

%% Description
% This function computes the predictive mean and variance at test input
% If SigX: consider random input, with covariance SigX. 
% Predictions computed using the exact equations. 
% The form of the covariance function is
%    C(x^p,x^q) = v1*exp( -0.5 * sum_{d=1..D} w_d * (x^p_d - x^q_d)^2 ) + v1
%    \delta_{pq}
% (computed using cov2.m)
% 
% Inputs: 
% loghteta .. optimized hyperparameters 
% covfunc .. dummy, used for (eventual) future compatibility  
% invQ .. inverse of the data covariance matrix
% input .. input part of the training data,  NxD matrix
% target .. output part of the training data (ie. target), Nx1 vector 
% muX       D by 1 test input
% SigX      covariance of the test input (OPTIONAL)
% lag .. the order of the model (number of used lagged outputs) 
% Outputs: 
% m  ..predicted mean 
% S2 .. predicted variance (noise free) 
%
% Based on the work of J. Quinonero-Candela and A. Girard. 

%% Examples
% demo_example_gp_simulation.m

%% See Also
% SIMUL00EXACT

[n, D] = size(input); % number of training cases and dimension of input space
[nn, D] = size(muX);  % number of test cases and dimension of input space
X=[-2*logtheta(1:end-2);2*logtheta(end-1);2*logtheta(end)]; % adapt hyperparameters to local format
expX = exp(X);        % exponentiate the hyperparameters once and for all

beta = invQ*target;

% Covariance between training and test inputs ...

a = zeros(n,nn);
for d = 1:D
  a = a + (repmat(input(:,d),1,nn)-repmat(muX(:,d)',n,1)).^2*expX(d);
end
a = expX(D+1)*exp(-0.5*a);

% Covariance between the test input and themselves 
b = expX(D+1);   

% Predictive mean  and variance (noise-free test input)
m = a'*beta;
S2 = b - sum(a.*(invQ*a),1)';

if  nargin > 6 % random test input
    
    L = sum(diag(SigX)~=0); % number of stochastic dimensions (non zero vars)
    if (L==0)
        return;
    end
    % little things to do before we start
    if nargin == 8 % in the case of control inputs
        rangeL = lag-L+1:lag;
        
    else
        rangeL = D-L+1:D;
        
    end
    rangeC = [1:rangeL(1)-1 rangeL(end)+1:D];
    
    SigX = SigX(rangeL,rangeL);
    muXL = muX(:,rangeL);
    muXC = muX(:,rangeC);
    inputL = input(:,rangeL);
    inputC = input(:,rangeC);
    
    invLL = diag(expX(rangeL));
    invLC = diag(expX(rangeC));
    invS = inv(SigX);
    invC = (invLL+invS);    
    invSmuX = invS*muXL';
    t1 = muXL*invSmuX;
    c = inv(invC)*(invLL*inputL'+repmat(invSmuX,1,n));
    t2 = sum(inputL.*(inputL*invLL),2);
    t3 = sum(c.*(invC*c),1)';
    I = (1/sqrt(det(invLL*SigX+eye(L))))*exp(-0.5*(t1+t2-t3));    
    CC = exp(-.5*(sum((inputC-repmat(muXC,n,1)).*((inputC-repmat(muXC,n,1))*invLC),2)));
    m = b.*(CC.*I)'*beta;
    
    invD = 2*invLL+invS;
    [kk1,kk2]=meshgrid(1:n);
    T1 = repmat(inputL,n,1)+inputL(reshape(kk1,1,n^2),:);
    invLT1 = T1*invLL;
    d = invD\(invLT1'+repmat(invSmuX,1,n^2));
    T3 = reshape(sum(d.*(invD*d),1),n,n);
    I2 = (1/sqrt(det(2*invLL*SigX+eye(L))))*exp(-0.5*(t1+repmat(t2,1,n)+repmat(t2',n,1)-T3));
    
    CCC = CC*CC';
    S2 = b - b^2*sum(sum((invQ-beta*beta').*(CCC.*I2))) - m^2;
   
end
