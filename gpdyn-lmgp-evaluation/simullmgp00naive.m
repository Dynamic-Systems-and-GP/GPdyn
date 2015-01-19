%% simullmgp00naive 
function [mu, s2] = simullmgp00naive(logtheta, covfunc, input, target, targetvariance,...
    derivinput, derivtarget, derivvariance, xt, lag)

%% Syntax
%  function [mu, s2] = simullmgp00naive(logtheta, covfunc, input, target, targetvariance,...
%   derivinput, derivtarget, derivvariance, xt, lag)

%% Description
% "Naive" (i.e. without propagation of variance) simulation of the GP
% model with the incorporated local models (LM). 
% Note: Currently it can be used only with Gaussian covariance function and
% with white noise model (sum of covSEard and covNoise). 
% Uses routine gpSD00. 
% 
% 
% Inputs: 
% loghteta .. optimized hyperparameters 
% covfunc .. specified covariance function, see help covFun for more info 
% input .. input part of the training data,  NxD matrix
% target .. output part of the training data (ie. target), Nx1 vector 
% targetvariance .. target variance, use NaN where not known 
% derivinput .. input part of the derivative training data, NEQxD matrix 
% derivtarget .. target derivatives, NEQxD matrix 
% derivvariance .. variances of the local model prameters, NEQxD matrix   
% xt .. input matrix for simulation, kxD vector, see
%   construct_simul_input.m for more info  
% lag .. the order of the model (number of used lagged outputs) 
% Outputs: 
% mu .. mean predicted output 
% s2 .. asociated variances 
% 
% Based on the work of R. Murray-Smith and A. Girard. 

%% Examples
% demo_example_lmgp_simulation.m

%% See Also
% GPSD00, SIMULLMGP00EXACT, SIMULLMGP00MCMC, SIMUL02NAIVE


fun_name = 'simullmgp00naive'; 

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



