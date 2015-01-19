function r = gmx_sample ( Mi, Sig, N )
% gmx_sample - creates samples of mixture components
%
%% Syntax
%  r = gmx_sample ( Mi, Sig, N )
%
%% Description
%  Creates samples of Gaussian mixture components
% 
% Input:
% * Mi    ... vector of mean values of mixture components 
% * Sig   ... vector of variances of mixture components
% * N     ... number of samples of mixture components
% 
% Output:
% * r     ... sample of the Gaussian mixture
% 
% See also: 
% simulGPnaive, gmx_sample
% 
% 
%% 
% * Written by J. Prikryl, November 2010

num_components = length(Mi); 
% First, create N uniform samples of mixture components.
%
% We will begin by determining the number of components in the mixture


% Then, uniformly sampling interval 1,num_components will provide us with
% the indices into Mi and Sig arrays. We are lucky as the weights of all
% mixture components in this case are uniform so we do not have to take
% them into account when sampling.
idx = randi ( num_components, [1,N] );

%
% Now sample the mixture components
%
r = randn(1,N).*Sig(idx)+Mi(idx);