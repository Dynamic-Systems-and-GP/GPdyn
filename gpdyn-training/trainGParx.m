 function [hyp, flogtheta, i] = trainGParx(hyp, inf, mean, cov, lik, input, target, minf)
% Function for optimization (training) of the GP model hyperparameters
% 
%% Syntax
%  [hyp, flogtheta, i] = trainGParx(hyp, inf, mean, cov, lik, input, target)
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
% * input  ... n by D matrix of training inputs
% * target ... column vector of length n of training targets
% * minf ... function handle of the minimization method to be used 
%                (optional, default=@minimize)
%
% Output:
% * hyp  ... optimized hyperparameters 
% * flogtheta ... minus log likelihood for the different runs (init. to 0)
% * i    ... number of iterations needed for the last optimization
%
% See Also:
% gp, minimize, trainlgmp, covFunctions, infMethods, likFunctions,
% meanFunctions 
% 
% Examples:
% demo_example_gp_training
%
%%

if(nargin < 8)
  minf = @minimize;
end

[n D] = size(input);

% logtheta0
% if(nargin<=3)        
%     logtheta = -1*rand(D+2,1); 
% else
%     logtheta = logtheta0; 
% end     

MIN_DIFF = 0.000001; %0.002; 
flogtheta = 0; 

[hyp, flogthetatmp,i] = feval(minf, hyp, @gp, -100, inf, mean, cov, lik, input, target);


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
    
    [hyp, flogthetatmp,i] = feval(minf,hyp, @gp, -100, inf, mean, cov, lik, input, target);
    
    if isempty(flogthetatmp) % no improvement: at minimum
        disp('oops')
        break
    end
    flogtheta = [flogtheta flogthetatmp(end)];
end



