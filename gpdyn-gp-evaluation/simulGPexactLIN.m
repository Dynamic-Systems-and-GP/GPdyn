function [m, s2, mu, sig2] = simulGPexactLIN(hyp, inf, mean, cov, lik, input, target, xt, lag)
% Simulation of the GP model, where the output variance is propagated using
% analytical approximation
%
%% Syntax
%  [mu, sig2, m, s2] = simulGPexactLIN(hyp, inf, mean, cov, lik, input, target, xt, lag);
%
%% Description
% See J.Kocijan, A. Girard, D.J. Leith, Incorporating linear local models 
% in Gaussian process model, Technical Report DP-8895, Institut Jozef Stefan, 
% Ljubljana, December 2003.
% Simulation of the GP model, where the output variance is propagated using
% analytical approximation. It can be used only with the linear 
% covariance function and with the white noise model (covLINard and likGauss). 
% Uses routine gpExactLINard. 
% 
% Input: 
% * hyp      ... the structure of optimized hyperparameters 
% * inf      ... the function specifying the inference method 
% * mean     ... the prior mean function
% * cov      ... the specified covariance function, see help covFun for more info 
% * lik      ... the likelihood function
% * input    ... the input part of the training data,  NxD matrix
% * target   ... the output part of the training data (ie. target), Nx1 vector 
% * xt       ... the input matrix for simulation, kxD vector, see
%                construct.m for more info 
% * lag      ... the the order of the model (number of used lagged outputs) 
% 
% Output: 
% * m     ... the predictive mean when propagating the uncertainty 
% * s2    ... the predictive variance when propagating the uncertainty (including noise variance)
% * mu    ... the predictive mean using "naive" approach (doesn't propagate the
%             uncertainty)
% * sig2  ... the predictive variance using "naive" approach   (including noise variance, as usual)
% 
% See also:
% covLINard, covNoise, gpExactLINard, simulGPExactSE
%
%% 
% * Written by J. Kocijan, 2010


fun_name = 'simulGPexactLIN'; 

[n, D] = size(input);
[nn, D] = size(xt);

%input validation
[ is_valid, hyp, inf, mean, cov, lik, msg ] = validate( hyp, inf, mean, cov, lik, D);

if ~isequal(cov,{@covLINard}) 
    error(strcat([fun_name,': function can only be called with the', ...
        ' covariance function ''covLINard'' '])); 
end

if ~isequal(lik,{@likGauss}) 
    error(strcat([fun_name,': function can only be called with the', ...
        ' likelihood function ''likGauss'', where hyp.lik parameter is log(sn)'])); 
end 

X=[-2*hyp.cov;2*hyp.lik];
x=input;
t=target;
expX = exp(X);

vy = expX(end);
Q = feval(cov{:},hyp.cov, x);
Q = Q + vy*eye(n);
invQ = inv(Q);


beta=invQ*t;
W=diag(exp(X(1:D)));
XW = x*W;
muXp = xt(1,:); % propagate uncertainty
SigX = zeros(D,D); 
muX = xt(1,:);  % naive approach

[mtmp, s2tmp] = gpExactLINard(hyp, inf, mean, cov, lik, invQ, input, target, muXp, SigX);
m(1) = mtmp; 
s2(1) = s2tmp;%

[mutmp, sig2tmp] = gpExactLINard(hyp, inf, mean, cov, lik, invQ, input, target, muX, SigX);
mu(1) = mutmp; 
sig2(1) = sig2tmp; 

for k=2:nn
    if (mod(k,50) == 0)
    disp(['simulGPexactLIN, step: ', int2str(k), '/', int2str(nn)]);
    end

    test=muXp;
    % For the NEXT prediction...
    
    % cross-cov terms
        
    gama=diag(exp(X(1:D)));
    ML_parameters=gama*x'*inv(x*gama*x'+eye(length(t),length(t))*exp(X(D+1)))*t; 
    covXY = (ML_parameters'*SigX)';
    % input covariance matrix and mean
    
    if nargin > 7 % control inputs to take into account
        SigX(1:lag-1,1:lag-1) = SigX(2:lag,2:lag);
        SigX(lag,lag) = s2(k-1);
        SigX(1:lag-1,lag) = covXY(2:lag);
        SigX(lag,1:lag-1) = covXY(2:lag)';
           
        muXp = [muXp(2:lag) m(k-1) xt(k,lag+1:end)];    
        muX = [muX(2:lag) mu(k-1) xt(k,lag+1:end)];
        [mtmp, s2tmp] = gpExactLINard(hyp, inf, mean, cov, lik, invQ, input, target, muXp, SigX);
    else
        SigX(1:D-1,1:D-1) = SigX(2:D,2:D);
        SigX(D,D) = s2(k-1);
        SigX(1:D-1,D) = covXY(2:D); 
        SigX(D,1:D-1) = covXY(2:D)'; 

        muXp = [muXp(2:D) m(k-1)];    
        muX = [muX(2:D) mu(k-1)];
        [mtmp, s2tmp] = gpExactLINard(hyp, inf, mean, cov, lik, invQ, input, target, muXp, SigX);
    end
        
    m(k) = mtmp; 
    s2(k) = s2tmp;
    [mutmp, sig2tmp] = gpExactLINard(hyp, inf, mean, cov, lik, invQ, input, target, muXp);
    mu(k) = mutmp; 
    sig2(k) = sig2tmp;    
    
    if s2tmp < 0
       keyboard
       
    end
    
end   

% transform into coloumn vectors 
m = m'; 
s2 = s2'; 
mu = mu';
sig2 = sig2'; 

% in the case of negative variance, but not to cover it!
for i=1:length(s2)
    if s2(i)<0
        s2(i)=0;
    end
    if sig2(i)<0
        sig2(i)=0;
    end
end


return; 

