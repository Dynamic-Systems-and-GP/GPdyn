function [y, s2] = predGPnaive(hyp, inf, mean, cov, lik, input, target, test, lag, kstep) 
% predGPNaive - 'Naive' (i.e. without propagation of variance) k-step-ahead prediction 
% of GP model.
%
%% Syntax
%   [y, s2] = predGPnaive(hyp, inf, mean, cov, lik, input, target, test, lag, kstep) 
%
%% Description
% GP-NARX model is used to predict a further step ahead by replacing the data at 
% present time instant with the data at one time instant before and using the mean 
% value of prediction from the previous prediction step instead of the measured 
% output value. This is then repeated for predifined number of steps, which is at 
% least 1.
% 
% Input: 
% * hyp      ... optimized hyperparameters (struct)
% * inf      ... the function specifying the inference method 
% * mean     ... the prior mean function
% * cov      ... the specified covariance function, see help covFun for more info 
% * lik      ... the likelihood function
% * input    ... the input part of the training data,  NxD matrix
% * target   ... the output part of the training data (ie. target), Nx1 vector 
% * test     ... the input matrix, kxD matrix, 
% * lag      ... the order of the model (number of used lagged outputs) 
% * kstep    ... the number of steps to be predicted
% 
% Output: 
% * y  ... the matrix of columns with mean outputs of predictions from 1 to kstep-ahead 
% * s2 ... the matrix of columns with associated variances 
%
% Examples:
% demo_example_gp_simulation
%
% See Also:
% gpx, simulGPnaive, simulGPmc
%

%% 
% * Written by J. Kocijan, 2009

[NN,D] = size(test);

% first step: k=1

    xt = test;  % naive approach

    [ytmp, s2tmp, fmu, fs, post] = gpx(hyp, inf, mean, cov, lik, input, target, xt);
    y = NaN(size(ytmp,1),kstep);                    % fill matrix with NaNs
    s2 = y;                                         % fill matrix with NaNs
    y(:,1) = ytmp;        % first step ahead prediction to the first column
    s2(:,1) = s2tmp;

% future steps ...   
for k=2:kstep
        
        % For the next prediction prepare the input vector ...
        if size(test,2)>lag
            Test = [xt(2:end,1:lag-1) y(k-1:end-1,k-1) xt(2:end,lag+1:D)];
            % Test =[delayed outputs, replaced last delayed output, delayed inputs] 
        else
            Test = [xt(2:end,2:D) y(k-1:end-1,k-1)];
        end
        
        xt=Test;          % ... and calculate the one-step-ahead prediction
        [ytmp, s2tmp, fmu, fs, post] = gpx(hyp, inf, mean, cov, lik, input, target, xt, post);
        y(k:end,k) = ytmp; 
        s2(k:end,k) = s2tmp;
        
        if s2tmp < 0
            keyboard
        end
        
    end   
end

