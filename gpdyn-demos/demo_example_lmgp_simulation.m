%% demo_example_lmgp_simulation 

%% Description
% Demo to present the simulation of the LMGP model of the chosen dynamic
% system. Three different simulations are presented: 
% (1) "naive" simulation (without propagation of uncertainty) 
% (2) "exact" simulation (anayltical propagation of uncertainty) 
% (3) simulation with numerical (MCMC) propagation of uncertainty. 
% See: 
%   A. Girard, Approximate Methods for Propagation of Uncertainty 
%   with Gaussian Process Models, PhD thesis, 2004. 
% for more info. 
% 
% Set flags accordingly. 
% 
% Note: 
% Currently it can be used only with Gaussian covariance function and
% with white noise model (sum of covSEard and covNoise). 

%% See Also
% EXAMPLE, DEMO_EXAMPLE_LMGP_DATA, DEMO_EXAMPLE_LMGP_TRAINING,
% SIMULLMGP00NAIVE, SIMULLMGP00EXACT, SIMULLMGP00MCMC, GPSD00, GPSD00RAN 

clear;
close all;

% validation data
load example_data

% trained LMGP model
load example_lmgp_trained

flag_naive_simulation = 1;
flag_exact_simulation = 1;
flag_mcmc_simulation = 1;

uvalid = example_valid_data.uvalid;
xvalid = example_valid_data.xvalid;
yvalid = example_valid_data.yvalid;

xt = [xvalid uvalid];
lag = 1;
t = [0:length(uvalid)-1]';


if(flag_naive_simulation)
    % naive simulation
    [ynaive, s2naive] = simullmgp00naive(logtheta, covfunc, input, target, targetvar,...
        inputDer, targetDer, derivevar, xt, lag); 

    plotgp(101,t,yvalid, ynaive, sqrt(s2naive));
    plotgpe(201,t,yvalid, ynaive, sqrt(s2naive));
end


if(flag_exact_simulation)
    % exact simulation
    [yexact, s2exact] = simullmgp00exact(logtheta, covfunc, input, target, targetvar,...
        inputDer, targetDer, derivevar, xt, lag); 


    plotgp(102,t,yvalid, yexact, sqrt(s2exact));
    plotgpe(202,t,yvalid, yexact, sqrt(s2exact));

end
 
if(flag_mcmc_simulation)
    % mcmc simulation
    Nsamples = 140;
    [ymcmc, s2mcmc, mcmcMM, mcmcVV] = simullmgp00mcmc(logtheta, covfunc, input,target,targetvar,...
        inputDer, targetDer, derivevar, xt, lag, Nsamples);

    plotgp(104,t,yvalid, ymcmc, sqrt(s2mcmc));
    plotgpe(204,t,yvalid, ymcmc, sqrt(s2mcmc));

    % testing probability density function in various steps
    steps_test_pdf = [1 135];
    mcmc_test_pdfs(mcmcMM,mcmcVV,[steps_test_pdf]);

end

return
