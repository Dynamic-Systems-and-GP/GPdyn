%% simul00mcmc
function [mu, s2, MU, SIG2, indexes_warning] = simul02mcmc(logtheta, covfunc, input, target, xt, lag, Nsamples)

%% Syntax
%  function [y, s2] = simul02mcmc(logtheta, covfunc, input, target, xt,lag) 

%% Description
% Simulation of the GP model, where the output variance is propagated using
% simple MCMC method, see A. Girard, Approximate Methods for Propagation of
% Uncertainty with Gaussian Process Models, PhD thesis, 2004. 
% Idea: at every time step the output of GP model is approximated with
% Nsamples samples, which are used as the future inputs of the GP model.
% Samples are re-used if necessary (ie. y(k-1) for y(k-2) if lag=2 etc.) 
% Uses routines gpr and mcmc_getsamplesgaussianmix. 
% 
% Inputs: 
% loghteta .. optimized hyperparameters 
% covfunc .. specified covariance function, see help covFun for more info 
% input .. input part of the training data,  NxD matrix
% target .. output part of the training data (ie. target), Nx1 vector 
% xt .. input matrix for simulation, kxD vector, see
%   construct_simul_input.m for more info 
% lag .. the order of the model (number of used lagged outputs) 
% Nsamples .. number of samples used in algorithm (ie. runs of simulation) 
% Outputs: 
% mu .. mean predicted output 
% s2 .. associated variances 
% MU .. matrix of all predicted means, kxNsamples
% SIG2 .. associated predicted variances 
% indexes_warning .. sampling problems, see mcmc_getsamplesgaussianmix.m
% 
% Written by K.Azman, 31.05.2005
% Based on the work of C.E. Rasmussen and A. Girard. 

%% Examples
% demo_example_gp_simulation.m

%% See Also
% GPR_SIMUL, MCMC_GETSAMPLESGAUSSIANMIX, SIMUL02NAIVE



Ndx = 800;
DSig = 3;
indexes_warning = [];

[N, D] = size(input);
PDF = zeros(Nsamples,lag);


% 1st step - input is a point
test = xt(1,:);
[mu(1), s2(1), alpha, L] = gpr_simul(logtheta, covfunc, input, target, test);

MU(1,:) = mu(1)*ones(1,Nsamples);
SIG2(1,:) = s2(1)*ones(1,Nsamples);

% 2nd step - input is Gaussian distribution
spacestart = mu(1)-DSig*sqrt(s2(1));
spaceend = mu(1)+DSig*sqrt(s2(1));
spacedx = (spaceend-spacestart)/Ndx;
inputspace = [spacestart:spacedx:spaceend]';

[dumm1,dumm2,dumm3,dumm4,pdf,fwarning] = mcmc_getsamplesgaussmix(inputspace,mu(1),s2(1),Nsamples);
PDF(:,lag) = pdf;
test = xt(2,:);
for ii=1:Nsamples
    test(lag) = pdf(ii);
    [MU(2,ii), SIG2(2,ii), alpha, L] = gpr_simul(logtheta, covfunc, input, target, test, alpha, L);
end
mu(2,1) = mean(MU(2,:));
s2(2,1) = mean(SIG2(2,:)) + mean((MU(2,:)-mu(2)).^2);

% steps from 3 on ...
for k=3:length(xt)

    if(mod(k,50)==0 | k==3)
        disp(['simul02mcmc, step: ',int2str(k),'/',int2str(length(xt))]);    
    end

    
    if(k>lag)
        % samples of previous distributions 
        for jj=1:lag-1
            PDF(:,jj) = PDF(:,jj+1);
        end

    else  % k <= lag
        % part after 'else' not yet tested for lag>=3 
        col0 = lag-(k-1);
        PDF(:,1:col0) = repmat(xt(k,1:col0),Nsamples,1);
        for jj=col0+1:lag-1
            PDF(:,jj) = PDF(:,jj+1);
        end
    end

    % generating current distribution
    spacestart = min(MU(k-1,:)-DSig*sqrt(SIG2(k-1,:)));
    spaceend = max(MU(k-1,:)+DSig*sqrt(SIG2(k-1,:)));
    spacedx = (spaceend-spacestart)/Ndx;
    inputspace = [spacestart:spacedx:spaceend]';

    [dumm1,dumm2,dumm3,dumm4,pdf, fwarning] = mcmc_getsamplesgaussmix(inputspace,MU(k-1,:),SIG2(k-1,:),Nsamples);
    PDF(:,lag) = pdf;

    if(fwarning==1)
        indexes_warning(end+1) = k;
    end

    % simulate for all cases
    for ii=1:Nsamples
        xx = PDF(ii,:);
        test =  [xx xt(k, lag+1:end)];
        [MU(k,ii), SIG2(k,ii), alpha, L] = gpr_simul(logtheta, covfunc, input, target, test, alpha, L);
    end

    % approximate output distribution with gauss - calculate m and v
    mu(k,1) = mean(MU(k,:));
    s2(k,1) = mean(SIG2(k,:)) + mean((MU(k,:)-mu(k)).^2);

end


return;





