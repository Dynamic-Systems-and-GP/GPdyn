%% simul02naive
function [y, s2] = simul02naive(logtheta, covfunc, input, target, xt, lag) 

%% Syntax
%  function [y, s2] = simul02naive(logtheta, covfunc, input, target, xt, lag) 

%% Description
% "Naive" (i.e. without propagation of variance) simulation of the GP
% model. 
% Uses routine gpr. 
% 
% Inputs: 
% loghteta .. optimized hyperparameters 
% covfunc .. specified covariance function, see help covFun for more info 
% input .. input part of the training data,  NxD matrix
% target .. output part of the training data (ie. target), Nx1 vector 
% xt .. input matrix for simulation, kxD vector, see
%   construct_simul_input.m for more info 
% lag .. the order of the model (number of used lagged outputs) 
% Outputs: 
% y .. mean predicted output 
% s2 .. associated variances 

%% Examples
% demo_example_gp_simulation.m

%% See Also
% GPR_SIMUL, SIMUL00EXACT, SIMUL02MCMC


[NN,D] = size(xt);

% first step: k=1
test = xt(1,:);
[y(1), s2(1), alpha, L] = gpr_simul(logtheta, covfunc, input, target, test);

% future steps ... 
for k=2:NN
    
    if (mod(k,100) == 0)
        disp(['simul02naive, step: ', int2str(k), '/', int2str(NN)]);
    end    
    
    if (k>lag)
        test = [y(k-lag:k-1) xt(k, lag+1:end)];
    elseif (k<=lag)
        test = [xt(k, 1:lag-k+1) y(1:k-1) xt(k, lag+1:end)];
    end


    [y(k), s2(k), alpha, L] = gpr_simul(logtheta, covfunc, input, target, test, alpha, L);
    
end


% transform into column vectors 
y = y'; 
s2 = s2'; 
