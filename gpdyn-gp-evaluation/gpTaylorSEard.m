function [mu, S2, deriv, S2deriv] = gpTaylorSEard(hyp, inf, mean, cov, lik, invQ, input, target, test, SigmaX)
% Computes the predictive meand and variance at the test input  and test inputs with 
% the Gaussian distribution for the squared exponential covariance function using the Taylor 
% approximation
% 
%% Syntax
% [m S2 deriv S2deriv] = gpTaylorSEard(hyp, inf, mean, cov, lik, invQ,
% input, target, test, SigmaX);
%
%% Description
% See A. Girard, Approximate Methods for Propagation of
% Uncertainty with Gaussian Process Models, PhD thesis, 2004.
% If input SigmaX exist, consider random input, i.e., input with the Gaussian distribution, 
% with covariance SigX. Predictions computed using the equations for 
% uncertainty propagation approximated with the Taylor linearisation. 
% It can be used only with the Gaussian covariance function 
% covSEard. Noise variance is added at the end. 
%
% Input:
% *  hyp      ... the structure of optimized hyperparameters 
% *  inf      ... the function specifying the inference method 
% *  mean     ... the prior mean function
% *  cov      ... the specified covariance function, see help covFun for more info 
% *  lik      ... the likelihood function
% *  input    ... the n by D matrix of training inputs
% *  target   ... the (column) vector (of size n) of targets
% *  test     ... the 1 by D test input
% *  SigmaX   ... is the D by D input covariance matrix (optional)
%
% Output:
% *  m       ... the (column) vector (of size nn) of prediced means
% *  S2      ... the (column) vector (of size nn) of predicted variances (including noise variance)
% *  deriv   ... the n by D matrix of mean partial derivatives
% *  S2deriv ... the n by D by D matrix of (co-)variances on the
%                partial derivatives (3 dimensional array (full covariance
%                per test point) 
%
% See Also:
% covSEard, covNoise, simulGPtaylorSE, gpExactSEard
%            
% Examples:
% demo_example_gp_simulation.m
%
%% Signature
% Written by K. Azman, 2007, based on the software of A. Girard. 

[n, D] = size(input);   % the number of training cases and dimension of input space
[nn, D] = size(test);       % the number of test cases and dimension of input space

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
expX = exp(X);              % exponentiate the hyperparameters once and for all

% Covariance between training and test inputs ...
a = zeros(n, nn);                                       % create and zero space
for d = 1:D
  a = a + expX(d)*(repmat(input(:,d),1,nn)-repmat(test(:,d)',n,1)).^2;
end
a = expX(D+1)*exp(-0.5*a);

% Covariance between the test input and themselves
b = expX(D+1) + expX(D+2); % signal variance plus noise variance

% Predicted mean  
invQt = invQ*target;
mu = a'*invQt; 

% Predicted variance
sig2 = b - sum(a.*(invQ*a),1)';
S2 = sig2;

if nargout > 2
    
    deriv = zeros(nn,D);
    S2deriv = zeros(D,D);
    for d = 1:D
        c = a.*(repmat(input(:,d),1,nn)-repmat(test(:,d)',n,1));
        deriv(1:nn,d) = expX(d)*c'*invQt; % derivative of mean
        ainvQc = a.*(invQ*c);
        for e = 1:D
            S2deriv(d,e) = -expX(d)*expX(e)* ...
                sum(ainvQc.*(repmat(input(:,e),1,nn)-repmat(test(:,e)',n,1)),1)';
        end
        S2deriv(d,d) = S2deriv(d,d) + expX(d)*expX(D+1);
    end

   if nargin == 10
    S2 = sig2 + deriv*SigmaX*deriv' + 0.5*sum(sum(SigmaX.*S2deriv));
   end
end
