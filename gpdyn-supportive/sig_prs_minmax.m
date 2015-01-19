%% sig_prs_minmax
function y = sig_prs_minmax(n, ksw, min, max)

%% Syntax
%  function y = sig_prs_minmax(n, ksw, min, max)

%% Description
%  Function, which generates "random" signal with length n and uniform
%  distribution of amplitudes between min and max. ksw is the minimum width
%  (duration in steps) of the signal at particular amplitude. 
% Inputs: 
% n .. length of the signal 
% ksw .. tact (minimum duration of signal at constant value) 
% min .. minimum amplitude 
% max .. maximum amplitude 
% Output: 
% y .. signal values

%% Examples
% demo_example_gp_data.m

k = ceil(n/ksw); 

rand('state',sum(100*clock)); 
y = rand(1,k); 

y = min + y*(max-min); 

y = repmat(y, ksw, 1); 
y = reshape(y, prod(size(y)),1); 



