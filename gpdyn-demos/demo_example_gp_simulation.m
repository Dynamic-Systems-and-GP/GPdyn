%% demo_example_gp_simulation 

%% Description
% This demo presents the simulation of GP model, describing dynamic
% system. Five different simulations are presented: 
% * "naive" simulation (without propagation of uncertainty) 
% * simulation with numerical (Monte Carlo) propagation of uncertainty. 
% * simulation using Taylor approximation for propagation of uncertainty. 
% * "exact" simulation (anayltical propagation of uncertainty) 
% * "exact" simulation (anayltical propagation of uncertainty) using linear
% covariance function
% 
% Based on work of: A. Girard, Approximate Methods for Propagation of
% Uncertainty with Gaussian Process Models, PhD thesis, 2004. 
 
%% See Also
% example, demo_example_present, demo_example_gp_data,
% demo_example_gp_training, demo_example_gp_norm, simulGPnaive,
% simulGPmc, simulGPtaylorSE, simulGPexactSE, simulGPexactLIN, construct

clear all;
close all;

% load data from file
load example_data 
load example_trained

% test data
y0 = 0; % initial value
lag = 1;

% construct input matrix for simulation (lag = 1):
% test = [y0  uvalid(1);
%          0  uvalid(2);
%
%          0  uvalid(end)];

test = construct(lag, uvalid, y0);  

t = [0:length(uvalid)-1]'; %time

%% Naive Simulation
[ynaive, s2naive] = simulGPnaive(hyp, inf, mean, cov, lik, input, target, test, lag);

f1=figure('Name', 'Naive Simulation');
plotgp(f1,t,yvalid, ynaive, sqrt(s2naive));

%% MC Simulation - Numerical Variance Propagation

Nsamples = 100;
[ymc, s2mc, MU, SIGMA2] = simulGPmc(hyp, inf, mean, cov, lik, input, target, test, lag,Nsamples);
f2=figure('Name', 'MC Simulation');
plotgp(f2,t,yvalid, ymc, sqrt(s2mc));

% testing pdfs in steps 1 and 14
desired_steps = [1 14]; 
mc_test_pdfs(MU,SIGMA2,desired_steps); 

%% Taylor Simulation - Variance Propagation using Taylor Approximation

[ytaylor, s2taylor] = ...
    simulGPtaylorSE(hyp, inf, mean, cov, lik, input, target, test, lag);

f4=figure('Name', 'Taylor Simulation');
plotgp(f4,t,yvalid, ytaylor, sqrt(s2taylor));

%% Exact Simulation - Analytical Variance Propagation

[yexactSE, s2exactSE] = ...
  simulGPexactSE(hyp, inf, mean, cov, lik, input, target, test, lag);

f3=figure('Name', 'Exact Simulation');
plotgp(f3,t,yvalid, yexactSE, sqrt(s2exactSE));

%% Exact Simulation - Analytical Variance Propagation (Using Linear Covariance Function)

[yexactLIN, s2exactLIN, ynaiveLIN, s2naiveLIN] = ...
    simulGPexactLIN(hyp_lin, inf, mean, @covLINard, lik, input, target, test, lag);

f4=figure('Name', 'Naive Simulation (Using Linear Covariance Function)');
plotgp(f4,t,yvalid, ynaiveLIN, sqrt(s2naiveLIN));
  
f5=figure('Name', 'Exact Simulation (Using Linear Covariance Function)');
plotgp(f5,t,yvalid, yexactLIN, sqrt(s2exactLIN));

%% Comparing results
% A function loss is available to calculate several frequently used
% perforamnce values. See the console output.

load example_trained
[yy, ss2] = simulGPexactSE(hyp, inf, mean, cov, lik, input, target, test, lag);

fprintf('\nNaive Simulation:\n');
loss(yvalid, ynaive, s2naive);                 
fprintf('\nMonte Carlo Simulation:\n');
loss(yvalid, ymc, s2mc);
fprintf('\nTaylor Simulation:\n');
loss(yvalid, ytaylor, s2taylor);
fprintf('\nExact Simulation:\n');
loss(yvalid, yexactSE, s2exactSE);  
fprintf('\nExact Simulation (Linear Covariance Function):\n');
loss(yvalid, yexactLIN, s2exactLIN);  



