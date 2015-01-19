function [ is_valid, hyp, inf, mean, cov, lik, msg ] = validate( hyp, inf, mean, cov, lik, D, throw_error)
% Function checks whether the hyperparameters and input functions match.
%
%% Syntax
% [ is_valid, hyp, inf, mean, cov, lik, msg ] = validate( hyp, inf, mean, cov, lik, D, throw_error)
%
%% Description
% Function checks if the the number of hiperparameters in the hyp struct
% matches the given covariance, likelihood and mean function. It also sets
% the default values for input functions if they are given as empty arrays.
% If hyp is also given as an empty array, its structure will not be checked,
% only the default input fuctions will be set it they are not given
% (@infExact, @meanZero and @likGauss) 
%
% This function is used for input validation in most of ruthines which
% perform regretion or simulation.
%
% Input:
% *  hyp  ... struct of hyperparameters
% *  inf  ... function specifying the inference method 
% *  cov  ... prior covariance function
% *  mean ... prior mean function
% *  lik  ... likelihood function
% *  D    ... size of the input space
% *  throw_error ... determines wether the error should be thrown (1) or
%                    returned as string. Default value is 1.
%
% Output:
% *  is_valid ... 1 if everything is ok, 0 on error
% *  hyp  ... struct of hyperparameters
% *  inf  ... function specifying the inference method 
% *  cov  ... prior covariance function (see below)
% *  mean ... prior mean function
% *  lik  ... likelihood function
% *  mgs  ... error message
%%
% * Written by Tomaž Šuštar

is_valid = 1;
msg = '';

% input validation
if nargin < 6
   error('Too few parameters are given.'); % 
end
if nargin < 7
    % deafaoul throw_error
    throw_error = 1;
end

if isempty(mean) %set default mean function
  mean=@meanZero;
end
if ischar(mean) || isa(mean, 'function_handle'), mean = {mean}; end  % make cell

if isempty(lik) %set default likelyhood function
  lik=@likGauss;
end
if ischar(lik) || isa(lik, 'function_handle'), lik = {lik}; end  % make cell

if isempty(inf) %set default inference method
  inf=@infExact;
end 
if ischar(inf) || isa(inf, 'function_handle'), inf = {inf}; end  % make cell

if isempty(cov) %covariance function can not be empty
  msg='Covariance function can not be empty';
  is_valid=0;
end
if ischar(cov) || isa(cov, 'function_handle'), cov = {cov}; end  % make cell

if ~isempty(hyp) % if hyp is empty do not chech hyperparameters
  if isstruct(hyp) % check format and size of hyp
    if ~isfield(hyp,'cov')
      hyp.cov=[];
    end
    if ~isfield(hyp,'mean')
      hyp.mean=[];
    end
    if ~isfield(hyp,'lik')
      hyp.lik=[];
    end

    hyp_cov_count = eval(eval_func(cov)); % number of covariance hyperparameters
    if length(hyp.cov) ~= hyp_cov_count
        is_valid = 0;
        msg = ['hyp.cov (size: ',num2str(length(hyp.cov)),') should have ', num2str(hyp_cov_count), ' hyperparameters.'] ;
    end

    hyp_lik_count = eval(eval_func(lik)); % number of likelihood hyperparameters
    if length(hyp.lik) ~= hyp_lik_count
        is_valid = 0;
        msg = [msg, ' hyp.lik (size: ',num2str(length(hyp.lik)),') should have ', num2str(hyp_lik_count), ' hyperparameters.'] ;
    end

    hyp_mean_count = eval(eval_func(mean)); % number of mean hyperparameters
    if length(hyp.mean) ~= hyp_mean_count
        is_valid = 0;
        msg = [msg,  ' hyp.mean (size: ',num2str(length(hyp.mean)),') should have ', num2str(hyp_mean_count), ' hyperparameters.'] ;
    end

  else
    is_valid = 0;
    msg = 'hyp should be a struct with fields hyp.cov hyp.lik and hyp.mean.';
  end
end
%throw error if needed
if throw_error && is_valid == 0
    error(msg);
end

end

