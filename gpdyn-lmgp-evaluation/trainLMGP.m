 function [hyp, nlZ, i] = trainLMGP(hyp0, inf, mean, cov, lik, input, target, inputDer, targetDer, deriveVar)
% trainLMGP - Function for optimization (training) of the LMGP model (GP model v local
% information) hyperparameters based on the training data via Maximum 
% Likelihood (ML).
%
%% Syntax
%  function [hyp, nlZ, i] = trainlmgp(hyp, inf, mean, cov, lik, input, target, inputDer, targetDer, deriveVar, hyp0)
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
% * hyp            ... a struct of hyperparameters
% x inf      	   ... the inference method 	  --> this is never used here
% x cov      	   ... prior covariance function  --> this is never used here
% x mean    	   ... prior mean function        --> this is never used here
% x lik      	   ... likelihood function        --> this is never used here
% * input          ... input part of the training data,  NxD matrix
% * target         ... output part of the training data (ie. target), Nx1 vector 
% * derivinput     ...input part of the derivative training data, NEQxD matrix  
% * derivtarget    ... target derivatives, NEQxD matrix 
% * derivevariance ... variances of the local model prameters, NEQxD matrix   
% * hyp0      ... intial hyperparameters (optional) 
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


[n D] = size(target);
[nd D] = size(targetDer);

  if ~isstruct(hyp0)
    hyp0.cov = -rand(D+1,1); 
    hyp0.lik = -rand(1,1); 
  end
nlZ = 0;

[hyp, nlZtmp, i] = minimize(hyp0, 'gpSD00', -200, inf, mean, cov, lik, input,target, NaN*ones(n,1),...
    inputDer,targetDer, deriveVar);

if isempty(nlZtmp)
    nlZ = nlZ;
else
    nlZ = [nlZ nlZtmp(end)];
end

while abs(nlZ(end) - nlZ(end-1)) > 0.001
	[hyp, nlZtmp, i] = minimize(hyp, 'gpSD00', -200, inf, mean, cov, lik, input,target,NaN*ones(n,1),...
		inputDer,targetDer, deriveVar);
	if isempty(nlZtmp) % no improvement: at minimum
	  break
	end
	nlZ = [nlZ nlZtmp(end)];
end
if exp(hyp.lik) < 1e-6
  hyp.lik =log(1e-6);
end




