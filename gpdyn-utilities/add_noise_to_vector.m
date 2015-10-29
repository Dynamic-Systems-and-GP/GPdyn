function [x_noise,x_no_noise] = add_noise_to_vector(x_no_noise, noise_std)
% Function adds the white Gaussian noise to vector given as a parameter. 
% 
%% Syntax 
% [x_noise,x_no_noise] = add_noise_to_vector(x_no_noise, noise_std)
% 
%% Description 
% It is used to get noisy data from simulation noise-free data. Function
% uses randn (internal Matlab command) to obtain sampled Gaussian
% distribution.  
% 
% Input:
% * x_no_noise ... the vector of noise-free data 
% * noise_std  ... the standard deviation of the white Gaussian noise 
% 
% Output:
% * x_noise    ... noisy data 
% * x_no_noise ... original data 
%
% 
% See also: 
% randn
% 
% Examples: demo_example_gp_data.m

noise = randn(size(x_no_noise));
noise = (noise - mean(noise))/std(noise)*noise_std;
x_noise = x_no_noise + noise; 

return; 