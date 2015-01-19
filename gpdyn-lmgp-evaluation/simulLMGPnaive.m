function [mu, s2] = simulLMGPnaive(logtheta, covfunc, input, target, targetvariance,...
    derivinput, derivtarget, derivvariance, xt, lag)
% simulLMGPnaive - 'Naive' (i.e. without propagation of variance) simulation of the GP
% model with the incorporated local models (LM).
%% Syntax
%  function [mu, s2] = simulLMGPnaive(logtheta, covfunc, input, target, targetvariance,...
%   derivinput, derivtarget, derivvariance, xt, lag)

%% Description
% The LMGP NARX model with the incorporated local models is used  
% to predict a further step ahead by replacing the data at 
% present time instant with the data at one time instant before and using the mean 
% value of prediction from the previous prediction step instead of the measured 
% output value. This is then repeated as many times as there are test samples.
% Currently it can be used only with Gaussian covariance function and
% with white noise model (sum of covSEard and covNoise). 
% Uses routine gpSD00. 
% Based on the work of R. Murray-Smith and A. Girard. 
% 
% Input: 
% * loghteta       ... optimized hyperparameters 
% * covfunc        ... specified covariance function, see help covFun for more info 
% * input          ... input part of the training data,  NxD matrix
% * target         ... output part of the training data (ie. target), Nx1 vector 
% * targetvariance ... target variance, use NaN where not known 
% * derivinput     ... input part of the derivative training data, NEQxD matrix 
% * derivtarget    ... target derivatives, NEQxD matrix 
% * derivvariance  ... variances of the local model prameters, NEQxD matrix   
% * xt             ... input matrix for simulation, kxD vector, see
%                      construct_simul_input.m for more info  
% * lag            ... the order of the model (number of used lagged outputs) 
%
% Output: 
% * mu             ... mean predicted output 
% * s2             ... asociated variances 
% 
% See Also
% gpSD00, simulLMGPmcmc, simulGPnaive
%
% Examples
% demo_example_LMGP_simulation.m

fun_name = 'simulLMGPnaive'; 

if ~(isequal(covfunc{1},'covSum') &  isequal(covfunc{2}{1},'covSEard') & ...
        isequal(covfunc{2}{2},'covNoise'))
    error(strcat([fun_name,': function can be called only with the sum', ...
        ' of covariance functions ''covSEard'' and ''covNoise'' '])); 
end 

% 1st point
test = xt(1,:); 
[mu(1), s2(1)] = gpSD00(logtheta, input, target, targetvariance, derivinput, derivtarget, derivvariance, test);
for k=2:size(xt,1)

    if(mod(k,50)==0)
        disp([fun_name,': ',int2str(k),'/',int2str(length(xt))]);
    end 

     
    if (k>lag)
        test = [mu(k-lag:k-1) xt(k, lag+1:end)];
    elseif (k<=lag)
        test = [xt(k, 1:lag-k+1) mu(1:k-1) xt(k, lag+1:end)];
    end
    


[mu(k), s2(k)] = gpSD00(logtheta, input, target, targetvariance,...
    derivinput, derivtarget, derivvariance, test);

end

% adding noise variance  
s2 = s2 + exp(2*logtheta(end));

mu=mu';
s2=s2';



