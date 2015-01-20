 function [logtheta, flogtheta, i] = trainLMGP(covfunc, likfunc, input, target, inputDer, targetDer, deriveVar, hyp0)
% trainLMGP - Function for optimization (training) of the LMGP model (GP model v local
% information) hyperparameters based on the training data via Maximum 
% Likelihood (ML).
%
%% Syntax
%  function [logtheta, flogtheta, i] = trainlmgp(covfunc, input, target, inputDer, targetDer, deriveVar, logtheta0)
%
%% Description
%  
% Function for optimization (training) of the LMGP model hyperparameters based on the training data via Maximum 
% Likelihood. It can be used only with Gaussian covariance function and
% with white noise model (sum of covSEard and covNoise). 
% Uses routines gpSD00 and minimize.
% Based on the work of C.E.Rasmussen and R. Murray-Smith.  
% 
% Input: 
% * covfunc        ... specified covariance function, see help covFun for more info 
% * input          ... input part of the training data,  NxD matrix
% * target         ... output part of the training data (ie. target), Nx1 vector 
% * derivinput     ...input part of the derivative training data, NEQxD matrix  
% * derivtarget    ... target derivatives, NEQxD matrix 
% * derivevariance ... variances of the local model prameters, NEQxD matrix   
% * logtheta0      ... intial values of hyperparameters (optional) 
%
% Output: 
% * loghteta       ... optimized hyperparameters 
% * flogtheta      ... minus log likelihood for the different runs (init. to 0)
% * i              ... number of iterations needed for the last optimization
%
% See Also
% minimize.m, gpSD00.m, trainGP.m
%
% Examples
% demo_example_LMGP_training.m


fun_name = 'trainLMGP'; 

if ~(isequal(covfunc{1},'covSEard'))
    error(strcat([fun_name,': function can be called only with covSEard'])); 
end 

if ~(isequal(likfunc{1},'likGauss'))
    error(strcat([fun_name,': function can be called only with likGauss'])); 
end 

[n D] = size(target);
[nd D] = size(targetDer);

  if nargin == 6
    hyp0 = struct('cov',-rand(D+2,1),'lik',-rand(),'mean',[]); 
  end
flogtheta = 0;

[logtheta, flogthetatmp, i] = minimize(hyp0, 'gpSD00', -200, input,target, NaN*ones(n,1),...
    inputDer,targetDer, deriveVar);

if isempty(flogthetatmp)
    flogtheta = flogtheta;
else
    flogtheta = [flogtheta flogthetatmp(end)];
end

while abs(flogtheta(end) - flogtheta(end-1)) > 0.001
[logtheta, flogthetatmp, i] = minimize(logtheta, 'gpSD00', -200, input,target,NaN*ones(n,1),...
    inputDer,targetDer, deriveVar);

if isempty(flogthetatmp) % no improvement: at minimum
  break
  end
  flogtheta = [flogtheta flogthetatmp(end)];
end
if exp(logtheta(end)) < 1e-6
  logtheta(end) =log(1e-6);
end




