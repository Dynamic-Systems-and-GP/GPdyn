function [ varargout ] = eval_func( func, varargin )
% Evaluates a function.
%
%% Syntax
% hyp_cov_count = eval(eval_func(cov));
%
%% Description
%   The robust method to evaluate covariance, mean and likelyhood functions in
%   the context of gpml toolbox. Input can be a function handle, string or
%   cell array.
% 
% Input:
% * func     ... the function handle, the function name (string) or cell array,
%             describing the function to be evaluated
% * varargin ... parameters of the input functon
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
% * Written by Tomaz Sustar, December 2011

if ischar(func) || isa(func, 'function_handle'), func = {func}; end  % make cell

varargout = {feval(func{:}, varargin)}; %evaluate function

end

