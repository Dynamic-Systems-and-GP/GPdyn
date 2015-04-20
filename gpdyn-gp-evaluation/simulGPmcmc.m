function [mu, s2, MU, SIG2] = simulGPmcmc(hyp, inf, meanfunc, cov, lik, input, target, test, lag, Nsamples)
% simulGPmcmc - Simulation of the dynamic GP model, where the output variance is
% propagated using simple MCMC method
%
%% Syntax
%  [mu, s2, MU, SIG2] = simulGPmcmc(hyp, inf, mean, cov, lik, input, target, test,
%  lag, Nsamples)
% 
%% Description
% Idea: at every time step the output of GP model is approximated with
% Nsamples samples, which are used as the future inputs of the GP model.
% Samples are re-used if necessary (ie. y(k-1) for y(k-2) if lag=2 etc.) 
% Uses routines gpx and gmx_sample. 
% 
% Input:
% * hyp      ... struct of optimized hyperparameters 
% * inf      ... function specifying the inference method 
% * meanfunc ... prior mean function
% * cov      ... specified covariance function, see help covFun for more info 
% * lik      ... likelihood function
% * input    ... input part of the training data,  NxD matrix
% * target   ... output part of the training data (ie. target), Nx1 vector 
% * test     ... input matrix for simulation, kxD vector, see
%                construct.m for more info  
% * lag      ... the order of the model (number of used lagged outputs) 
% * Nsamples ... number of samples used in algorithm (ie. runs of simulation) 
% 
% Output:
% * mu    ... mean predicted output 
% * s2    ... associated variances (with noise variances)
% * MU    ... matrix of all predicted means, kxNsamples
% * SIG2  ... associated predicted variances 
% 
% See also: 
% gpx, gmx_sample, simulGPnaive
% 
% Examples: 
% demo_example_gp_simulation
% 
%% 
% * Written by J. Prikryl, November 2010
% * Based on the work of C.E. Rasmussen, A. Girard, K. Azman. 
%

% meanfunc ... 'mean' is used as a matlab core function in this file


Ndx = 800;
DSig = 3;
num_iters = length(test);
sum_time  = 0;

[N, D] = size(input);
PDF = zeros(Nsamples,lag);

% Preallocate mu and s2
mu = zeros ( num_iters, 1 );
s2 = zeros ( num_iters, 1 );

% 1st step - input is a point
test_ = test(1,:);
[mu(1), s2(1), post] = gpx(hyp, inf, meanfunc, cov, lik, input, target, test_);

MU(1,:) = mu(1)*ones(1,Nsamples);
SIG2(1,:) = s2(1)*ones(1,Nsamples);

pdf = gmx_sample(mu(1),s2(1),Nsamples);
PDF(:,lag) = pdf;
% instead of a cycle, create one monstrous test_ matrix
test_ = repmat( test(2,:), Nsamples, 1 );
test_(:,lag) = pdf;
[MU(2,:), SIG2(2,:), post] = gpx(hyp, inf, meanfunc, cov, lik, input, target, test_, post);
mu(2,1) = mean(MU(2,:));
s2(2,1) = mean(SIG2(2,:)) + mean((MU(2,:)-mu(2)).^2);

% steps from 3 on ...
for k=3:num_iters

    if(mod(k,50)==0 || k==3)
        disp(['simulGPmcmc, step: ',int2str(k),'/',int2str(length(test))]);    
    end

    
    if(k>lag)
        % samples of previous distributions 
        for jj=1:lag-1
            PDF(:,jj) = PDF(:,jj+1);
        end

    else  % k <= lag
        % part after 'else' not yet tested for lag>=3 
        col0 = lag-(k-1);
        PDF(:,1:col0) = repmat(test(k,1:col0),Nsamples,1);
        for jj=col0+1:lag-1
            PDF(:,jj) = PDF(:,jj+1);
        end
    end

    pdf = gmx_sample(MU(k-1,:),SIG2(k-1,:),Nsamples);
    PDF(:,lag) = pdf;

    % simulate for all cases, again using the matrix version
    t0 = tic;
    test_ = repmat ( test(k,:), Nsamples, 1 );
    test_(:,1:lag) = PDF;
    [MU(k,:), SIG2(k,:), post] = gpx(hyp, inf, meanfunc, cov, lik, input, target, test_, post);
    calltime = toc(t0);
%     fprintf ( 'gpr_simul() ....... %f sec\n\n', calltime ); % uncomment
%     if you wish 
    sum_time = sum_time + calltime;
    
    % approximate output distribution with gauss - calculate m and v
    mu(k,1) = mean(MU(k,:));
    s2(k,1) = mean(SIG2(k,:)) + mean((MU(k,:)-mu(k)).^2);

end

%     fprintf ( 'average time needed for gpr_simul() .... %f sec\n', sum_time/(num_iters-2) );
%     uncomment if you wish
return;





