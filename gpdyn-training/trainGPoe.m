function [hyp, flogtheta, i] = trainGPoe(hyp, inf, mean, cov, lik, input, target, simf, lag, Nsamples)
% Function for the optimisation (training) of GP model hyperparameters
%
%% Syntax
% [hyp, flogtheta, i] = trainGPoe(hyp, inf, mean, cov, lik, input, target,
%                                 simf, lag, Nsamples);
%
%% Description
% Function for the optimisation (training) of GP model hyperparameters
% based on the training data via maximum marginal likelihood. 
% Uses routines gp and minimize.
% Based on the work of C.E.Rasmussen. 
% 
% Input:
% * hyp      ... the structure of initial hyperparameters
% * inf      ... the function specifying the inference method 
% * cov      ... the prior covariance function (see below)
% * mean     ... the prior mean function
% * lik      ... the likelihood function
% * input    ... the input part of the training data,  NxD matrix
% * target   ... the output part of the training data (ie. target), Nx1 vector
% * simf     ... the function handle of simulation function (e.g. @simulGPmc)
% * lag      ... the order of the model (number of used lagged outputs)
% * Nsamples ... the number of samples for MCMC simulation (optional)
% * minf     ... the function handle of the minimization method to be used 
%                (optional, default=@minimize)
%
% Output: 
% * hyp       ... optimized hyperparameters 
% * flogtheta ... the minus log likelihood for the different runs (init. to 0)
% * i         ... the number of iterations needed for the last optimization
%
% Examples:
% demo_example_gp_training.m
%
% See Also:
% gp, minimize, covFunctions, trainlgmp
%

if (nargin < 10)
  Nsamples = 100; % default value of Nsamples for MC simulation
end
if(nargin < 11)
  minf = @minimize;
end

if (nargin < 9)
  error('Too few parameters are given.');
end

y = feval(simf, hyp, inf, mean, cov, lik, input, target, input, lag, Nsamples);

inputSim = input(lag:end,:);
for i = 1 : lag
  inputSim(:,i) = y(lag+1-i:end+1-i);
end

MIN_DIFF = 0.002; 
flogtheta = 0; 

[hyp, flogthetatmp, i] = feval(minf, hyp, @gp, -100, inf, mean, cov, lik, inputSim, target);

if isempty(flogthetatmp)
  error('Minimization failed - matrix close to singular.');
else
  flogtheta = [flogtheta flogthetatmp(end)];
end

while (abs(flogtheta(end) - flogtheta(end-1))>MIN_DIFF)         
  disp(' '); 
  disp(strcat(['delta flogtheta: ', num2str(abs(flogtheta(end) - flogtheta(end-1)))])); 
  disp(' ')
  
  [logtheta, flogthetatmp,i] = feval(minf, hyp, @gp, -100, inf, mean, cov, lik, input, target);
    
  if isempty(flogthetatmp) % no improvement: at minimum
    disp('oops');
    break
  end
  flogtheta = [flogtheta flogthetatmp(end)];
end
