function [ varargout ] = eval_func( func, varargin )
% Evaluates a function.
%
%% Syntax
% hyp_cov_count = eval(eval_func(cov));
%
%% Description
%   Robust method to evaluate covariance, mean and likelyhood functions in
%   the context of gpml toolbox. Input can be a function handle, string or
%   cell array.
% 
% Input:
% * func  ... Fuction handle, funtion name (string) or cell array,
%             desribing a function to be evaluated
% * varargin ... Parameters of the input functon
% 
% Output:
% This function returns whatever is returned by the evaluated function
%
% See Also:
% validate
%
% Examples:
% validate
%
%%
% * Written by Tomaž Šuštar, December 2011

if ischar(func) || isa(func, 'function_handle'), func = {func}; end  % make cell

varargout = {feval(func{:}, varargin)}; %evaluate function

end

