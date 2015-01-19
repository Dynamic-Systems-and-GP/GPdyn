%% traingp
 function [logtheta, flogtheta, i] = traingp(covfunc, input, target, logtheta0)

%% Syntax
%  function [logtheta, flogtheta, i] = traingp(covfunc, input, target, logtheta0)

%% Description
% Function for optimization (training) of the GP model hyperparameters
% based on the training data via Maximum Likelihood (ML). 
% Uses routines gpr and minimize.
% Based on the work of C.E.Rasmussen. 
% 
% Inputs: 
% covfunc .. specified covariance function, see help covFun for more info 
% input .. input part of the training data,  NxD matrix
% target .. output part of the training data (ie. target), Nx1 vector 
% logtheta0 .. intial values of hyperparameters (optional) 
% Outputs: 
% logtheta .. optimized hyperparameters 
% flogtheta .. minus log likelihood for the different runs (init. to 0)
% i .. number of iterations needed for the last optimization

%% Examples
% demo_example_gp_training.m

%% See Also
% GPR, MINIMIZE, COVFUN, TRAINLMGP


[n D] = size(input);

% logtheta0
if(nargin<=3)        
    logtheta = -1*rand(D+2,1); 
else
    logtheta = logtheta0; 
end     

MIN_DIFF = 0.002; 
flogtheta = 0; 

[logtheta, flogthetatmp,i] = minimize(logtheta, 'gpr', -100, covfunc, input, target);


if isempty(flogthetatmp)
    flogtheta = flogtheta;
    error('minimization failed - matrix close to singular');
else
    flogtheta = [flogtheta flogthetatmp(end)];
end


while (abs(flogtheta(end) - flogtheta(end-1))>MIN_DIFF)     

    disp(' '); 
    disp(strcat(['delta flogtheta: ', num2str(abs(flogtheta(end) - flogtheta(end-1)))])); 
    disp(' ')
    
    [logtheta, flogthetatmp,i] = minimize(logtheta, 'gpr', -100, covfunc, input, target);
    
    if isempty(flogthetatmp) % no improvement: at minimum
        disp('oops')
        break
    end
    flogtheta = [flogtheta flogthetatmp(end)];
end



