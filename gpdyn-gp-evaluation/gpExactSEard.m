function [m, S2] = gpExactSEard(hyp, inf, mean, cov, lik, invQ, input, target, muX, SigX, lag)
% Computes the predictive mean and variance at teh test input with
% Gaussian distribution for the squared exponential covariance function.
%
%% Syntax
%  [m, S2] = gpExactSEard(hyp, inf, mean, cov, lik, invQ, input, target, muX, SigX, lag);
%
%% Description
% See A. Girard, Approximate Methods for Propagation of
% Uncertainty with Gaussian Process Models, PhD thesis, 2004.
% If input SigX exist, consider random input, i.e., input with Gaussian distribution, 
% with covariance SigX. Predictions computed using the equations for 'exact' 
% approximation of uncertainty propagation. It can be used only 
% with the Gaussian covariance function covSEard. Noise variance is added at the end.
% 
% Input: 
% * hyp      ... optimized hyperparameters 
% * inf      ... the function specifying the inference method 
% * mean     ... the prior mean function
% * cov      ... the specified covariance function, see help covFun for more info 
% * lik      ... the likelihood function
% * invQ     ... the inverse of the data covariance matrix
% * input    ... the input part of the training data,  NxD matrix
% * target   ... the output part of the training data (ie. target), Nx1 vector 
% * muX      ... the D by 1 test input
% * SigX     ... the covariance of the test input (OPTIONAL)
% * lag      ... the order of the model (number of used lagged outputs) 
%
% Output: 
% * m  ... the predicted mean 
% * S2 ... the predicted variance  (including noise variance)
%
% See Also:
% covSEard, simulGPexactSE, gpTaylorSEard
%
% Examples:
% demo_example_gp_simulation.m
%
%% Signature
% Written by K. Azman, 2007, based on the software of J. Quinonero-Candela 
% and A. Girard.

[n, D] = size(input); % the number of training cases and dimension of input space
[nn, D] = size(muX);  % the number of test cases and dimension of input space

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


X=[-2*hyp.cov(1:end-1);2*hyp.cov(end);2*hyp.lik]; % adapt hyperparameters to local format
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

% Predictive mean  and variance (test input including noise variance)
m = a'*beta;
S2 = b - sum(a.*(invQ*a),1)'  + expX(D+2);

if  nargin > 9 % random test input
    
    L = sum(diag(SigX)~=0); % number of stochastic dimensions (non zero vars)
    if (L==0)
        return;
    end
    % little things to do before we start
    if nargin == 11 % in the case of control inputs
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
    S2 = b - b^2*sum(sum((invQ-beta*beta').*(CCC.*I2))) - m^2 + expX(D+2); % including noise variance
   
end
