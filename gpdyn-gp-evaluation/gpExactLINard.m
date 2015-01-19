function [m, S2] = gpExactLINard(hyp, inf, mean, cov, lik, invQ, input, target, test, SigmaX)
% This function computes the predictive mean and variance at test input
% with Gaussian distribution for linear covariance function.
%
%% Syntax
%  [m, S2] = gpExactLINard(hyp, inf, mean, cov, lik, invQ, input, target, test, SigmaX);
%
%% Description
% See J.Kocijan, A. Girard, D.J. Leith, Incorporating linear local models 
% in Gaussian process model, Technical Report DP-8895, Institut Jožef Stefan, 
% Ljubljana, December 2003.
% This function computes the predictive mean and variance at test input
% If input SigmaX exist consider random input, with covariance SigmaX. 
% Predictions computed using the equations for 'exact' approximation 
% of uncertainty propagation. It can be used only 
% with linear covariance function covLINard. covNoise is not accounted
% for.
% 
% Inputs: 
% * hyp      ... struct of optimized hyperparameters 
% * inf      ... function specifying the inference method 
% * mean     ... prior mean function
% * cov      ... specified covariance function, see help covFun for more info 
% * lik      ... likelihood function
% * invQ     ... inverse of the data covariance matrix
% * input    ... input part of the training data,  NxD matrix
% * target   ... output part of the training data (ie. target), Nx1 vector 
% * test     ... D by 1 test input
% * SigmaX   ... covariance of the test input (OPTIONAL)
% 
% Outputs: 
% * m  ... predicted mean 
% * S2 ... predicted variance (noise free) 
%
%
% See also:
% covLINard, covNoise, simulGPexactLIN
% 
%% 
% Written by J. Kocijan, 2010

[n, D] = size(input); % number of training cases and dimension of input space
[nn, D] = size(test);  % number of test cases and dimension of input space

%input validation
% [ is_valid, hyp, inf, mean, cov, lik, msg ] = validate( hyp, inf, mean, cov, lik, D);
% 
% if ~isequal(cov,{@covLINard}) 
%     error(strcat([fun_name,': function can only be called with the', ...
%         ' covariance function ''covLINard'' '])); 
% end
% 
% if ~isequal(lik,{@likGauss}) 
%     error(strcat([fun_name,': function can only be called with the', ...
%         ' likelihood function ''likGauss'', where hyp.lik parameter is log(sn)'])); 
% end 

X=[-2*hyp.cov;2*hyp.lik]; % adapt hyperparameters to local format

expX = exp(X);        % exponentiate the hyperparameters once and for all

% ... then we compute the covariance between training and test inputs ...

a = zeros(n, nn);                                       % create and zero space
for d = 1:D                                               % linear contribution
  a = a + expX(d)*input(:,d)*test(:,d)';
end


% ... and covariance between the test input and themselves 

b = test.^2*expX(1:D);


% Now, write back mean prediction and variance

m = a'*invQ*target;
% Predicted variance
sig2 = b - sum(a.*(invQ*a),1)';
S2 = sig2;
if nargin > 9
    W=diag(expX(1:D));
    alpha=W-(input*W)'*invQ*(input*W);
    beta=invQ*target;
    gamma=(input*W)'*beta*beta'*(input*W);
    S2=test*alpha*test'+trace(alpha*SigmaX)+trace(gamma*SigmaX);
end
