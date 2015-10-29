 function [hyp, flogtheta, i] = trainGParx(hyp, inf, mean, cov, lik, input, target, minf)
% Function for the optimisation (training) of GP model hyperparameters
% 
%% Syntax
%  [hyp, flogtheta, i] = trainGParx(hyp, inf, mean, cov, lik, input, target)
% 
%% Description
% Function for the optimisation (training) of GP model hyperparameters
% based on the training data via maximum marginal likelihood. 
% Uses routines gp and minimize.
% Based on the work of C.E.Rasmussen. 
% 
% Input:
% * hyp    ... the structure of initial hyperparameters
% * inf    ... the function specifying the inference method 
% * cov    ... the prior covariance function (see below)
% * mean   ... the prior mean function
% * lik    ... the likelihood function
% * input  ... the n by D matrix of training inputs
% * target ... the column vector of length n of training targets
% * minf   ... the function handle of the minimization method to be used 
%                (optional, default=@minimize)
%
% Output:
% * hyp       ... optimized hyperparameters 
% * flogtheta ... the minus log likelihood for the different runs (init. to 0)
% * i         ... the number of iterations needed for the last optimization
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



