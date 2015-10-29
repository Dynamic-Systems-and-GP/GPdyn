function [y, s2, fmu, fs2] = simulGPnaive(hyp, inf, mean, cov, lik, input, target, test, lag) 
% simulGPnaive - 'Naive' (i.e. without propagation of variance) simulation of the GP
% ARX and AR model.
%
%% Syntax
%   [y, s2] = simulGPnaive(hyp, inf, mean, cov, lik, input, target, test, lag) 
%
%% Description
% GP-NARX model is used to predict a further step ahead by replacing the data at 
% present time instant with the data at one time instant before and using the mean 
% value of prediction from the previous prediction step instead of the measured 
% output value. This is then repeated as many times as there are test samples.
%
% Input:
% * hyp    ... the column vector of hyperparameters
% * inf    ... the function specifying the inference method
% * cov    ... the prior covariance function (see below)
% * mean   ... the prior mean function
% * lik    ... the likelihood function
% * input  ... the input part of the training data,  NxD matrix
% * target ... the output part of the training data (ie. target), Nx1 vector 
% * test   ... the input matrix, kxD matrix, see construct.m for more info 
% * lag    ... the order of the model (number of used lagged outputs) 
% 
% Output:
% * y    ... the mean predicted output 
% * s2   ... associated variances
% * fmu  ... the mean predicted output without noise
% * fs2  ... associated variances without noise
%
% See also: 
% gpx.m, simulGPmc, construct.m
%
% Examples: 
% demo_example_gp_simulation.m
%
%% 
% Written by K. Azman in J. Kocijan, 2007

[NN,D] = size(test); % D - input space dimension, NN - number or test samples

% allocate memory
y = zeros(1,D);
s2 = y;
fmu = y;
fs2 = y;

% first step: k=1
xt = test(1,:);
[y(1), s2(1) ,fmu(1), fs2(1), post] = gpx(hyp, inf, mean, cov, lik, input, target, xt);

% future steps ... 
for k=2:NN
    
    if (mod(k,100) == 0)  % remark on every 100th step
        disp(['simulGPnaive, step: ', int2str(k), '/', int2str(NN)]);
    end
        % For the next prediction prepare the input vector ...
        if D > lag
            if (k>lag)
                xt = [y(k-lag:k-1) test(k, lag+1:end)];
            elseif (k<=lag)
                xt = [test(k, 1:lag-k+1) y(1:k-1) test(k, lag+1:end)];
            end
        else                                            % if D <= lag
            if (k>lag)
                xt = y(k-lag:k-1);
            elseif (k<=lag)
                xt = [test(k, 1:lag-k+1) y(1:k-1)];
            end
        end
        % ... and calculate the one-step-ahead prediction
    [y(k), s2(k) ,fmu(k), fs2(k), post] = gpx(hyp, inf, mean, cov, lik, input, target, xt, post);
    
end


% transform into column vectors 
y = y'; 
s2 = s2'; 
fmu = fmu';
fs2 = fs2';
