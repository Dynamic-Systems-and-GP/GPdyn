function [mu, s2] = simulLMGPnaive(hyp, inf, mean, cov, lik, input, target, targetvariance,...
    derivinput, derivtarget, derivvariance, xt, lag)
% simulLMGPnaive - 'Naive' (i.e. without propagation of variance) simulation of the GP
% model with the incorporated local models (LM).
%% Syntax
%  function [mu, s2] = simulLMGPnaive(hyp, inf, mean, cov, lik, input, target, targetvariance,...
%   derivinput, derivtarget, derivvariance, xt, lag)

%% Description
% The LMGP NARX model with the incorporated local models is used  
% to predict a further step ahead by replacing the data at 
% present time instant with the data at one time instant before and using the mean 
% value of prediction from the previous prediction step instead of the measured 
% output value. This is then repeated as many times as there are test samples.
% Currently it can be used only with covSEard covariance function and
% with likGauss likelihood. 
% Uses routine gpSD00.
% Based on the work of R. Murray-Smith and A. Girard. 
% 
% Input: 
% * hyp            ... a struct of hyperparameters
% x inf      	   ... the inference method 	  --> this is never used here
% x cov      	   ... prior covariance function  --> this is never used here
% x mean    	   ... prior mean function        --> this is never used here
% x lik      	   ... likelihood function        --> this is never used here
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
% gpSD00, simulLMGPmc, simulGPnaive
%
% Examples
% demo_example_LMGP_simulation.m
%
% Changelog:
%
% 16.2.2015, Martin Stepancic:
%		 	-changed the function interface as gpml > 3.0
%			-removed the addition of autocovariance - this is now
%			 already included in gpSD00.m
%
fun_name = 'simulLMGPnaive'; 


% 1st point
test = xt(1,:); 
[mu(1), s2(1)] = gpSD00(hyp, inf, mean, cov, lik, input, target, targetvariance, derivinput, derivtarget, derivvariance, test);
for k=2:size(xt,1)

    if(mod(k,50)==0)
        disp([fun_name,': ',int2str(k),'/',int2str(length(xt))]);
    end 

     
    if (k>lag)
        test = [mu(k-lag:k-1) xt(k, lag+1:end)];
    elseif (k<=lag)
        test = [xt(k, 1:lag-k+1) mu(1:k-1) xt(k, lag+1:end)];
    end
    


[mu(k), s2(k)] = gpSD00(hyp, inf, mean, cov, lik, input, target, targetvariance,...
    derivinput, derivtarget, derivvariance, test);

end

mu=mu';
s2=s2';











