function [hyp, flogtheta0] = gp_initial(bounds, inf, mean, cov, lik, input, target, islog, npop)
% Function for finding initial values of hyperparameters with the random
% search.
%
%% Syntax
%  [hyp, flogtheta0] = gp_initial(bounds, inf, mean, cov, lik, input,
%  target, islog, npop)
%
%% Description 
% Returns best set of n random sets of hyperparameter values. As score it 
% uses a log marginal likelihood.
% 
% Input:
% * bounds ... bounds values of hyperaparameters, vector [min, max] if all
%              hyperparameters have the same bounds, otherwise vector of
%              two hyperparameter structs [hyp_min, hyp_max], if bouns
%              equals [] default bounds are used [-8,7]
% * inf    ... the function specifying the inference method 
% * mean   ... the prior mean function
% * cov    ... the prior covariance function
% * lik    ... the likelihood function
% * input  ... the input part of the training data,  NxD matrix
% * target ... the output part of the training data (ie. target), Nx1 vector

% * islog  ... indicates if random values should be in log-scale [optional,
%             default = 1]
% * npop   ... the number of population (sets of hyperparameters) [optional]
%
% Output: 
% * hyp        ... the structure icluding the initial hyperparameters 
% * flogtheta0 ... the log marginal likelihood of initial hyperparameters
%
% See Also:
% likelyhood, gp.m, minimize, covFunctions, likFunctions, meanFunctions,
% infMethod
%
% Examples:
% demo_example_gp_training.m
%
%% 
% * Written by Dejan Petelin (2010-05-03)
% * Modified by Tomaz Sustar, March 2012


if nargin < 7 % input validation
  error('Too few parameters are given.'); % 
end


% log-scale and number of population
if nargin < 8
  islog = 1;
end
% number of population
if nargin < 9    
    npop = 1000; % default 1000 init values
end

% dimensions
    [n D] = size(input);

% default values
  [ is_valid, hyp, inf, mean, cov, lik] = validate( [], inf, mean, cov, lik, D);
  
% dimensions
    
    hyp_cov_count = eval(eval_func(cov)); % number of covariance hyperparameters
    hyp_lik_count = eval(eval_func(lik)); % number of likelihood hyperparameters
    hyp_mean_count = eval(eval_func(mean)); % number of mean hyperparameters
    hyp_count = hyp_cov_count + hyp_lik_count + hyp_mean_count;
   
    % create hyperparameter struct
    hyp.cov=zeros(hyp_cov_count,1);
    hyp.lik=zeros(hyp_lik_count,1);
    hyp.mean=zeros(hyp_mean_count,1);
    
    
% bounds

% default bounds, log-scale and number of population
if isempty(bounds)
  bounds = repmat([-8 7],hyp_count,1); % default bound from -8 to 7
elseif length(bounds) == 2
    if isnumeric(bounds)
        bounds = repmat(bounds,hyp_count,1);
    elseif validate( bounds(1), inf, mean, cov, lik, D) && ...
            validate( bounds(2), inf, mean, cov, lik, D)

        bounds = [unwrap(bounds(1)), unwrap(bounds(2))];
        
    else
        error('Given bounds are not in a corret format!');
    end    
end

if (~islog && (min(min(bounds)) < 0))
  error('If random values should be in linear scale, then bound values must be equal or greater than zero!');
end

% random values
diff = bounds(:,2)-bounds(:,1);
logthetas = repmat(bounds(:,1),1,npop)+rand(hyp_count,npop).*repmat(diff,1,npop);

% if linear scale, values must be logarithmized
% for example: in GPR hyperparameters are given as log("linear value") or
%              as "logarithmic value", see gpr_demo for details
if (~islog)
  logthetas = log(logthetas);
end

% back to structs
%allocate
hyps = repmat(hyp, 1, size(logthetas,2));


for i = 1:size(logthetas,2)
    hyps(i)=rewrap(hyp, logthetas(:,i));
end

% negative marginal likelihoods
ml = likelihood(hyps, inf, mean, cov, lik, input, target);

% best
[flogtheta0 mi] = min(ml);
hyp = hyps(mi);

end


