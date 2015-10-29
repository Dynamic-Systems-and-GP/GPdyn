function [ml] = likelihood(hyp, inf, mean, cov, lik, input, target)
% Calculates the negative log marginal likelihood. 
%
%% Syntax
%  ml = likelihood(hyp, inf, mean, cov, lik, input, target)
%
%% Description
% Function for calculating the negative log marginal likelihood of given
% hyperparameters, covariance function and data. If more sets of
% hyperparameters are given it outputs a vector of calculated negative log
% marginal likelihoods.
%
% Based on the work of C.E.Rasmussen. 
% 
% Input: 
% * hyp    ... the hyperparameter struct(s), row vector of structs
% * inf    ... the function specifying the inference method 
% * mean   ... the prior mean function
% * cov    ... the prior covariance function
% * lik    ... the likelihood function
% * input  ... the input part of the training data,  NxD matrix
% * target ... the output part of the training data (ie. target), Nx1 vector
%
% Output:
% * ml     ... the negative log marginal likelihood(s), vector
%
% See Also:
% gpx, minimize, covFunction, gp_initial
%
% Examples:
% gp_initial.m
%%

if nargin < 7 % input validation
  error('Too few parameters are given.'); % 
end

ml=zeros(size(hyp)); % allocate

for i = 1 : size(hyp,2)
  try
    
    mlt = gp(hyp(i), inf, mean, cov, lik, input, target);
  catch
    mlt = +Inf;
  end

  % stability
  if (~isfinite(mlt) || isnan(mlt) || ~isreal(mlt))
    mlt = +Inf;
  end
  ml(i) = mlt;
end