% demo_example_lmgp_simulation - demonstration script file for LMGP model identification - simulation part. 
%
%% Description
% Demo to present the simulation of the LMGP model of the chosen dynamic
% system. Three different simulations are presented: 
% (1) 'naive' simulation (without propagation of uncertainty) 
% (2) simulation with numerical (Monte Carlo) propagation of uncertainty. 
% See: 
%   A. Girard, Approximate Methods for Propagation of Uncertainty 
%   with Gaussian Process Models, PhD thesis, 2004. 
% for more info. 
% 
% Set flags accordingly. 
%  
% Currently it can be used only with Gaussian covariance function and
% with white noise model (sum of covSEard and covNoise). 
%
% See Also
% example.m, demo_example_LMGP_data.m, demo_example_LMGP_training.m,
% simulLMGPnaive.m, simulLMGPmc, gpSD00 

clear;
close all;

% validation data
mat_data=load('example_data_lmgp_oeq.mat')

% trained LMGP model
load example_lmgp_trained


flag_naive_simulation = 1;
flag_mc_simulation = 1;

uvalid = mat_data.valid_data.u;
xvalid = mat_data.valid_data.x;
yvalid = mat_data.valid_data.y;

xt = [xvalid uvalid];
lag = 1;
t = [0:length(uvalid)-1]';


if(flag_naive_simulation)
    % naive simulation
    [ynaive, s2naive] = simulLMGPnaive(hyp, inffunc, meanfunc, covfunc, likfunc , ...
        input, target, targetvar, inputDer, targetDer, derivevar, xt, lag); 

    plotgp(101,t,yvalid, ynaive, sqrt(s2naive));
    plotgpe(201,t,yvalid, ynaive, sqrt(s2naive));
end


if(flag_mc_simulation)
    % mc simulation
    Nsamples = 40;
    [ymc, s2mc, mcMM, mcVV] = simulLMGPmc(hyp, inffunc, meanfunc, ...
        covfunc, likfunc , input,target,targetvar,inputDer, targetDer, derivevar, ...
        xt, lag, Nsamples);

    plotgp(104,t,yvalid, ymc, sqrt(s2mc));
    plotgpe(204,t,yvalid, ymc, sqrt(s2mc));

    % testing probability density function in various steps
    steps_test_pdf = [1 135];
    mc_test_pdfs(mcMM,mcVV,[steps_test_pdf]);

end

return
