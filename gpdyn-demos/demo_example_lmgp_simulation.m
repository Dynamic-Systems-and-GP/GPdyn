% demo_example_lmgp_simulation - demonstration script file for LMGP model identification - simulation part. 
%
%% Description
% Demo to present the simulation of LMGP model of the chosen dynamic
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
% Currently it can be used only with the Gaussian covariance function and
% with the white noise model (sum of covSEard and covNoise). 
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
    
    f1=figure('Name', 'Naive simulation');
    plotgp(f1,t,yvalid, ynaive, sqrt(s2naive));
    f2=figure('Name', 'Naive simulation');
    plotgpe(f2,t,yvalid, ynaive, sqrt(s2naive));
end


if(flag_mc_simulation)
    % MC simulation
    Nsamples = 100;
    [ymc, s2mc, mcMM, mcVV] = simulLMGPmc(hyp, inffunc, meanfunc, ...
        covfunc, likfunc , input,target,targetvar,inputDer, targetDer, derivevar, ...
        xt, lag, Nsamples);
    
    f3=figure('Name', 'Monte Carlo simulation');
    plotgp(f3,t,yvalid, ymc, sqrt(s2mc));
    f4=figure('Name', 'Monte Carlo simulation');
    plotgpe(f4,t,yvalid, ymc, sqrt(s2mc));

    % testing probability density function in various steps
    steps_test_pdf = [1 135];
    mc_test_pdfs(mcMM,mcVV,[steps_test_pdf]);

end

return
