%% trainlmgp
 function [logtheta, flogtheta, i] = trainlmgp(covfunc, input, target, inputDer, targetDer, deriveVar, logtheta0)

%% Syntax
%  function [logtheta, flogtheta, i] = trainlmgp(covfunc, input, target, inputDer, targetDer, deriveVar, logtheta0)

%% Description
% Function for optimization (training) of the LMGP model (GP model v local
% information) hyperparameters based on the training data via Maximum 
% Likelihood (ML). 
% Note: currently it can be used only with Gaussian covariance function and
% with white noise model (sum of covSEard and covNoise). 
% Uses routines gpSD00 and minimize.
% % Based on the work of C.E.Rasmussen and R. Murray-Smith.  
% 
% Inputs: 
% covfunc .. specified covariance function, see help covFun for more info 
% input .. input part of the training data,  NxD matrix
% target .. output part of the training data (ie. target), Nx1 vector 
% derivinput ..input part of the derivative training data, NEQxD matrix  
% derivtarget .. target derivatives, NEQxD matrix 
% derivevariance .. variances of the local model prameters, NEQxD matrix   
% logtheta0 .. intial values of hyperparameters (optional) 
% Outputs: 
% loghteta .. optimized hyperparameters 
% flogtheta .. minus log likelihood for the different runs (init. to 0)
% i .. number of iterations needed for the last optimization

%% Examples
% demo_example_lmgp_training.m

%% See Also
% MINIMIZE, GPSD00, TRAINGP

fun_name = 'trainlmgp'; 

if ~(isequal(covfunc{1},'covSum') &  isequal(covfunc{2}{1},'covSEard') & ...
        isequal(covfunc{2}{2},'covNoise'))
    error(strcat([fun_name,': function can be called only with the sum', ...
        ' of covariance functions ''covSEard'' and ''covNoise'' '])); 
end 

[n D] = size(target);
[nd D] = size(targetDer);

logtheta0 = -rand(D+2,1); 
flogtheta = 0;

[logtheta, flogthetatmp, i] = minimize(logtheta0, 'gpSD00', -200, input,target, NaN*ones(n,1),...
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




