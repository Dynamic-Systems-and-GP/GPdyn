function out = sig_prs_minmax(n, ksw, min, max)
%
%% Syntax
%  function out = sig_prs_minmax(n, ksw, min, max)
%
%% Description
%  Function, which generates pseudo random signal with length n and uniform
%  distribution of signal magnitudes between min and max. ksw is the minimum 
%  duration in samples of the signal at particular magnitude.
%
% Input: 
% * n   ... length of the signal 
% * ksw ... minimum number of samples of signal at constant value 
% * min ... minimum magnitude 
% * max ... maximum magnitude 
%
% Output: 
% * out ... PRS signal values
%
% See also:
% sig_prbs
%
%% Examples
% demo_example_gp_data.m
%
%% 
% * Written by K. Azman, 2007  


k = ceil(n/ksw); 

rand('state',sum(100*clock)); 
out = rand(1,k); 

out = min + out*(max-min); 

out = repmat(out, ksw, 1); 
out = reshape(out, prod(size(out)),1); 



