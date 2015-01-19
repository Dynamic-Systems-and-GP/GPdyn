%% simullmgp00exact
 function [mu, s2] = simullmgp00exact(logtheta, covfunc, input, target, targetvar,...
    derivinput, derivtarget, derivvariance, xt, lag)

%% Syntax
%  [mu, s2] = simullmgp00exact(logtheta, covfunc, input, target, targetvar,...
%    derivinput, derivtarget, derivvariance, xt, lag)

%% Description
% Simulation of the GP model with incorporated local models (LMGP models),
% where the output variance is propagated using analytical approximation,
% see J. Kocijan, A. Girard, and D. J. Leith. Incorporating linear local
% models in Gaussian process model. Technical Report DP-8895,
% Institut Jožef Stefan, Ljubljana, December 2003.
% A. Girard, Approximate Methods for Propagation of Uncertainty with
% Gaussian Process Models, PhD thesis, 2004. 
% Notes: 
% Currently it can be used only with Gaussian covariance function and
% with white noise model (sum of covSEard and covNoise). 
%
% Uses routine gpSD00ran. 
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
% Based on the work of J. Kocijan, A. Girard, R. Murray-Smith. 

%% Examples
% demo_example_lmgp_simulation.m

%% See Also
% GPSD00RAN, SIMULLMGP00NAIVE, SIMULLMGP00MCMC, SIMUL00EXACT



fun_name = 'simullmgp00exact'; 

if ~(isequal(covfunc{1},'covSum') &  isequal(covfunc{2}{1},'covSEard') & ...
        isequal(covfunc{2}{2},'covNoise'))
    error(strcat([fun_name,': function can be called only with the sum', ...
        ' of covariance functions ''covSEard'' and ''covNoise'' '])); 
end 


[n, D]      = size(input); 
fullinput  = [input; repmat(derivinput,D,1)];
fulltarget = [target; derivtarget(:)];
[N, D]      = size(fullinput); 
[nD, D] = size(derivinput); 
SigmaX = 0.0*eye(2*lag,2*lag);

% 1st point
test = xt(1,:); 

[mu(1), s2(1), invQ] = gpSD00ran(logtheta, input, target, targetvar,...
    derivinput, derivtarget, derivvariance, test, SigmaX);



W=inv(diag(exp(logtheta(1:D))));
beta=invQ*fulltarget;
for k=2:length(xt)
    
    if(mod(k,50)==0)
        disp([fun_name,': ',int2str(k),'/',int2str(length(xt))]);
    end
    
    
    % For the NEXT prediction...
    
    % cross-cov terms
    [nn, D] = size(test);
    a = zeros(N, 1);    % compute the covariance between training and test cases
    
    V=inv(W+SigmaX);
    for d=1:D
        for e = 1:D
            a = a + (fullinput(:,d)-repmat(test(:,e)',N,1)).^2*V(d,e); 
        end
    end
    
    
    a = exp(2*logtheta(D+1))*det(eye(D)+inv(W)*SigmaX)^(-0.5)*exp(-0.5*a);
    C=W-W*inv(W+SigmaX)*W;
    c=fullinput*(eye(D)-inv(W+SigmaX)*W)+repmat(test,N,1)*inv(inv(W)*SigmaX+eye(D));   
    ad=[];
    ac=[];
    for d=1:D
        ac(1:n,d)=a(1:n).*c(1:n,d);
        ad=[ad;
            -exp(logtheta(d)).*a(n+(d-1)*nD+1:d*nD+n,1)*...
                (c(n+(d-1)*nD+1:d*nD+n,:)'*fullinput(n+(d-1)*nD+1:d*nD+n,d)-c(n+(d-1)*nD+1:d*nD+n,:)'*c((d-1)*nD+1:d*nD,d)-C(:,d))'];
    end
    ac=[ac;ad];
    covXY=(beta'*ac)'-mu(k-1)*test';
    
    
    % input covariance matrix and mean
    
    SigmaX(1:lag-1,1:lag-1) = SigmaX(2:lag,2:lag);
    SigmaX(lag,lag) = s2(k-1);
    SigmaX(1:lag-1,lag) = covXY(2:lag);
    SigmaX(lag,1:lag-1) = covXY(2:lag)';
    
         
    if (k>lag)
        test = [mu(k-lag:k-1) xt(k, lag+1:end)];
    elseif (k<=lag)
        test = [xt(k, 1:lag-k+1) mu(1:k-1) xt(k, lag+1:end)];
    end
    
    [mu(k), s2(k), invQ] = gpSD00ran(logtheta, input, target, targetvar,...
        derivinput, derivtarget, derivvariance, test, SigmaX);
    
end   

% adding noise variance  
s2 = s2 + exp(2*logtheta(end));

mu = mu'; 
s2 = s2'; 


