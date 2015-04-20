function [hyp, flogtheta, i] = trainGPoe(hyp, inf, mean, cov, lik, input, target, simf, lag, Nsamples)
% Function for optimization (training) of the GP model hyperparameters
%
%% Syntax
% [hyp, flogtheta, i] = trainGPoe(hyp, inf, mean, cov, lik, input, target,
%                                 simf, lag, Nsamples);
%
%% Description
% Function for optimization (training) of the GP model hyperparameters
% based on the training data via Maximum Likelihood (ML). 
% Uses routines gp and minimize.
% Based on the work of C.E.Rasmussen. 
% 
% Input:
% * hyp  ... struct of initial hyperparameters
% * inf  ... function specifying the inference method 
% * cov  ... prior covariance function (see below)
% * mean ... prior mean function
% * lik  ... likelihood function
% * input    ... input part of the training data,  NxD matrix
% * target   ... output part of the training data (ie. target), Nx1 vector
% * simf     ... function handle of simulation function (e.g. @simulGPmc)
% * lag      ... the order of the model (number of used lagged outputs)
% * Nsamples ... number of samples for MCMC simulation (optional)
% * minf     ... function handle of the minimization method to be used 
%                (optional, default=@minimize)
%
% Outputs: 
% * hyp       ... optimized hyperparameters 
% * flogtheta ... minus log likelihood for the different runs (init. to 0)
% * i         ... number of iterations needed for the last optimization
%
% Examples:
% demo_example_gp_training.m
%
% See Also:
% gp, minimize, covFunctions, trainlgmp
%
%% Signature

if (nargin < 10)
  Nsamples = 100; % default value of Nsamples for MCMC simulation
end
if(nargin < 11)
  minf = @minimize;
end

% if (nargin < 5)
%   [n D] = size(input);
%   logtheta0 = -rand(D+2, 1);
% end
if (nargin < 9)
  error('Too few parameters are given.');
end

y = feval(simf, hyp, inf, mean, cov, lik, input, target, input, lag, Nsamples);

inputSim = input(lag:end,:);
for i = 1 : lag
  inputSim(:,i) = y(lag+1-i:end+1-i);
  %inputSim(i,:) = y(lag+1-i:end+1-i);
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
