
function [mu, s2, MM, VV] = simulLMGPmc(hyp, inffcn, meanfcn, covfcn, likfcn, input, target, targetvariance,...
    derivinput, derivtarget, derivvariance, xt, lag, Nsamples)
% simulLMGPmc - Simulation of the dynamic GP model with incorporated local models (LMGP models),
% where the output variance is propagated using numerical Monte Carlo approximation
%
%% Syntax
%  [mu, s2, MM, VV] = simulLMGPmc(hyp, inf, mean, cov, lik, input, target, targetvariance,...
%     derivinput, derivtarget, derivvariance, xt, lag, Nsamples)
%
%% Description
% Idea: at every time step the output of GP model is approximated with
% the Gaussian distribution, from which we sample one value. We repeat the
% procedure Nsamples-times. 
% Nsamples samples, which are used as the future inputs of the GP model.
% Samples are re-used if necessary (ie. y(k-1) for y(k-2) if lag=2 etc.) 
% Currently it can be used only with Gaussian covariance function and
% with white noise model (sum of covSEard and covNoise) due to the gpSD00. 
% Uses routine gpSD00. 
% see K. Azman. Identifikacija dinami�nih sistemov z Gaussovimi procesi. PhD
% thesis, Univerza v Ljubljani, Ljubljana, 2007. (in Slovene). 
%  
% Input: 
% * hyp            ... a struct of hyperparameters
%   inf      	   ... the inference method 	  --> this is never used here
%   cov      	   ... prior covariance function  --> this is never used here
%   mean    	   ... prior mean function        --> this is never used here
%   lik      	   ... likelihood function        --> this is never used here
% * input          ... input part of the training data,  NxD matrix
% * target         ... output part of the training data (ie. target), Nx1 vector 
% * targetvariance ... target variance, use NaN where not known 
% * derivinput     ... input part of the derivative training data, NEQxD matrix 
% * derivtarget    ... target derivatives, NEQxD matrix 
% * derivvariance  ... variances of the local model prameters, NEQxD matrix   
% * xt             ... input matrix for simulation, kxD vector, see
%                      construct_ARXsimul_input.m for more info 
% * lag            ... the order of the model (number of used lagged outputs) 
% * Nsamples       ... number of samples used in algorithm (ie. runs of simulation)
%
% Output: 
% * mu             ... mean predicted output 
% * s2             ... associated variances (including noise variance)
% * MM             ... matrix of all predicted means, kxNsamples
% * VV             ... associated predicted variances 
%
% See Also
% gpSD00.m, simulLMGPnaive.m, simulGPmc.m
%
% Examples
% demo_example_lmgp_simulation.m
%
%%
% * Written by K.Azman, 31.05.2005
% * Based on the work of C.E. Rasmussen and A. Girard. 
%
%
% Changelog:
%
% 16.2.2015, Martin Stepancic:
%		 	-changed the function interface as gpml > 3.0
%			-removed the addition of autocovariance - this is now
%			 already included in gpSD00.m. Beware: till now, the underlying
%			 p.d.f. for sampling did not contain the autocovariance term!
%





fun_name = 'simulLMGPmc'; 



[n, D]      = size(input);
fullinput  = [input; repmat(derivinput,D,1)];
fulltarget = [target; derivtarget(:)];
[N, D]      = size(fullinput);
[nD, D] = size(derivinput);


MM = zeros(Nsamples,size(xt,1));
VV = zeros(Nsamples,size(xt,1));

for jjj=1:Nsamples

    if(mod(jjj,10)==0)
        disp([fun_name, ', run:  ',int2str(jjj),'/',int2str(Nsamples)]);
    end

    % 1st point - input is "point" 
    test = xt(1,:);
    [mu(1), s2(1)] = gpSD00(hyp, inffcn, meanfcn, covfcn, likfcn, input, target, ...
			targetvariance, derivinput, derivtarget, derivvariance, test);

    for k=2:length(xt)

        test = [xt(k, lag+1:end)];

        % For the NEXT prediction...
        % assumed normal distribution, for more accurate procedure see
        % simulGPmc 

        % random sample from assumed distribution 
        ysampled = mu(k-1) + randn(1)*sqrt(s2(k-1));
        % generate input for the GP model 
        if (k>lag)
            test = [mu(k-lag:k-2) ysampled xt(k, lag+1:end)];
        elseif (k<=lag)
            test = [xt(k, 1:lag-k+1) mu(1:k-2) ysampled xt(k, lag+1:end)];
        end

        [mu(k), s2(k)] = gpSD00(hyp, inffcn, meanfcn, covfcn, likfcn, input, target, ...
			targetvariance, derivinput, derivtarget, derivvariance, test);
    end

    MM(jjj,:) = mu;
    VV(jjj,:) = s2;

end
% individual realisations saved in matrices MM and VV 
% approximate all output distributions with Gaussian distribution 
mu = mean(MM);
s2 = mean(VV) + mean((MM-repmat(mu,Nsamples,1)).^2);

mu = mu';
s2 = s2';


MM = MM'; 
VV = VV'; 

return 














