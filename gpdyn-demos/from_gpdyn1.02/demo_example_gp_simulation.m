%% demo_example_gp_simulation 

%% Description
% Demo to present the simulation of the GP model, describing dynamic
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

%% See Also
% EXAMPLE, DEMO_EXAMPLE_PRESENT, DEMO_EXAMPLE_GP_DATA, 
% DEMO_EXAMPLE_GP_TRAINING, SIMUL02NAIVE, SIMUL00EXACT, SIMUL02MCMC, GPR_SIMUL 

clear;
close all;

flag_naive_simulation = 1;
flag_exact_simulation = 1;
flag_mcmc_simulation = 1;

if (abs(flag_naive_simulation)+abs(flag_exact_simulation)+abs(flag_mcmc_simulation) == 0)
    disp('demo_example_simulation: at least one simulation should be chosen');
end

load example_data
load example_trained

uvalid = example_valid_data.uvalid;
xvalid = example_valid_data.xvalid;
yvalid = example_valid_data.yvalid;

% test input 
xt = [xvalid uvalid];
lag = 1;
t = [0:length(uvalid)-1]';



% naive simulation
if(flag_naive_simulation)
    [ynaive, s2naive] = simul02naive(logtheta,covfunc,input,target,xt,lag);
    plotgp(101,t,yvalid, ynaive, sqrt(s2naive));
    plotgpe(201,t,yvalid, ynaive, sqrt(s2naive));
end




% exact simulation - analytical variance propagation
if(flag_exact_simulation)
    [ynaive_exact, s2naive_exact, yexact, s2exact] = ...
        simul00exact(logtheta,covfunc,input,target,xt,lag);

    plotgp(102,t,yvalid, yexact, sqrt(s2exact));
    plotgp(103,t,yvalid, ynaive_exact, sqrt(s2naive_exact));

    plotgpe(202,t,yvalid, yexact, sqrt(s2exact));
    plotgpe(203,t,yvalid, ynaive_exact, sqrt(s2naive_exact));
end

% mcmc simulation - numerical variance propagation
if(flag_mcmc_simulation)
    Nsamples = 100;
    [ymcmc, s2mcmc, MU, SIGMA2] = simul02mcmc(logtheta,covfunc,input,target,xt,lag,Nsamples);
    
    plotgp(104,t,yvalid, ymcmc, sqrt(s2mcmc));
    plotgpe(204,t,yvalid, ymcmc, sqrt(s2mcmc));
    
    % testing pdfs in steps 1 and 14
    desired_steps = [1 14]; 
    mcmc_test_pdfs(MU,SIGMA2,desired_steps); 

end


return



