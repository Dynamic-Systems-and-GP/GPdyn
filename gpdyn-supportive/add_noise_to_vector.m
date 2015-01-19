%% add noise to vector 
function [x_noise,x_no_noise] = add_noise_to_vector(x_no_noise, noise_std)

%% Syntax
% function [x_noise,x_no_noise] = add_noise_to_vector(x_no_noise, noise_std)

%% Description
% Function adds white Gaussian noise to vector given as a parameter. 
% It is used to get noisy data from simulation noise-free data. 
% Inputs: 
% x_no_noise .. vector of noise-free data 
% noise_std .. standard deviation of white Gaussian noise 
% Outputs: 
% x_noise .. noisy data 
% x_no_noise .. original data 
% 
% Function uses randn (internal Matlab command) to obtain sampled Gaussian 
% distribution. 

%% Examples
% demo_example_gp_data.m

%% See Also
%
% RANDN

noise = randn(size(x_no_noise));
noise = (noise - mean(noise))/std(noise)*noise_std;
x_noise = x_no_noise + noise; 

return; 