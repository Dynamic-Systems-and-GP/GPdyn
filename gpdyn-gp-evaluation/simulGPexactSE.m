function [m, s2, mu, sig2] = simulGPexactSE(hyp, inf, mean, cov, lik, input, target, test, lag)
% Simulation of the GP model, where the output variance is propagated using
% the 'exact' analytical approximation
%
%% Syntax
%  [m, s2, mu, sig2] = simulGPexactSE(hyp, inf, mean, cov, lik, input, target, test, lag)
%
%% Description
% See A. Girard, Approximate Methods for Propagation of
% Uncertainty with Gaussian Process Models, PhD thesis, 2004. 
% Simulation of GP model, where the output variance is propagated using
% the 'exact' analytical approximation. It can be used only with the Gaussian 
% covariance function and the Gaussian likelihood function with the white noise
% model (covSEard and likGauss).  
% Uses routine gpExactSEard. 
% 
% Input: 
% *  hyp      ... the structure of optimized hyperparameters 
% *  inf      ... the function specifying the inference method 
% *  mean     ... the prior mean function
% *  cov      ... the specified covariance function, see help covFun for more info 
% *  lik      ... the likelihood function
% *  input    ... the input part of the training data,  NxD matrix
% *  target   ... the output part of the training data (ie. target), Nx1 vector 
% *  test     ... the input matrix for simulation, kxD vector, see
%                 construct.m for more info 
% *  lag      ... the order of the model (number of used lagged outputs) 
% 
% Output: 
% *  m     ... the predictive mean when propagating the uncertainty 
% *  s2    ... the predictive variance when propagating the uncertainty    (including noise variance)
% *  mu    ... the predictive mean using 'naive' approach (doesn't propagate
%              the uncertainty)
% *  sig2  ... the predictive variance using 'naive' approach   (including noise variance, as usual)
% 
% See also:
% gpExactSEard, simulGPtaylorSE, simulGPnaive
%
% Examples:
% demo_example_gp_simulation.m
%

%% 
% Written by K. Azman, 2007, based on the software of J. Quinonero-Candela 
% and A. Girard. 

fun_name = 'simulGPexactSE'; 


[n, D] = size(input);
[nn, D] = size(test);

% input validation
[ is_valid, hyp, inf, mean, cov, lik, msg ] = validate( hyp, inf, mean, cov, lik, D);

if ~isequal(cov,{@covSEard}) 
    error(strcat([fun_name,': function can only be called with the', ...
        ' covariance function ''covSEard'' '])); 
end

if ~isequal(lik,{@likGauss}) 
    error(strcat([fun_name,': function can only be called with the', ...
        ' likelihood function ''likGauss'', where hyp.lik parameter is log(sn)'])); 
end 


%define the lagged signal each regressor (system output and control input) represents:
if nargin==9
	switch class(lag)
	case 'double'
		assert(numel(lag)==1);
		assert(lag>0);
		maxlag=lag;
		ylags=(lag:-1:1);
		ulags=(lag:-1:1);
	case 'cell'
		maxlag=max([lag{:}]);
		ylags=lag{1};
		ulags=lag{2};	
	end
else
	ylags=(D:-1:1);
	ulags=[];
end
lags=[ylags,ulags+maxlag];


X=[-2*hyp.cov(1:end-1);2*hyp.cov(end);2*hyp.lik]; % adapt hyperparameters to local format
expX = exp(X);

vy = expX(end);


% training covariance matrix
Q = zeros(n,n);
for d = 1:D
  Q = Q + (repmat(input(:,d),1,n)-repmat(input(:,d)',n,1)).^2*expX(d);
end
Q = expX(D+1)*exp(-0.5*Q);

Q = Q + vy*eye(n);  % data cov. matrix: add noise
invQ = inv(Q);
alphaQ = invQ*target; %beta really
b = expX(D+1); % v_1 really...

post=infExact(hyp,mean,cov,lik,input,target);
post.Q=Q;
post.invQ=invQ;
post.alphaQ=alphaQ;


% 1st point
muXp = [test(1,1:maxlag) test(1,maxlag+1:2*maxlag)]; % mean values of regressors for prediction with propagation
SigX = zeros(2*maxlag,2*maxlag);% zero covariance matrix of regressors (assuming noise-free regressors are given)
muX = muXp;  % mean values of regressors for naive prediction (without propagation)

[mtmp, s2tmp, vtmp] = predictExactSEard(hyp,input,muXp(lags),SigX(lags,lags),post);
m(1) = mtmp; 
s2(1) = s2tmp; 

[mutmp, sig2tmp] = gpx(hyp, inf, mean,cov,lik, input, target, muX(lags),post);
mu(1) = mutmp; 

sig2(1) = sig2tmp; 

for k=2:nn
    if (mod(k,50) == 0)
        disp(strcat([fun_name, ', step: ', int2str(k), '/', int2str(nn)]));
    end

    SigX(2:maxlag,2:maxlag) = SigX(1:maxlag-1,1:maxlag-1);
    SigX(1,1) = s2(k-1);
    ioc=SigX(lags,lags)*vtmp;
    SigX(lags(ylags~=1),1) = ioc(1:length(ylags)-1);
    SigX(1,lags(ylags~=1)) = ioc(1:length(ylags)-1)';
    
    muXp = [m(k-1) muXp(1:maxlag-1) test(k,maxlag+1:2*maxlag)];    
    muX = [mu(k-1) muX(1:maxlag-1) test(k,maxlag+1:2*maxlag)];
    
    [mtmp, s2tmp, vtmp] = predictExactSEard(hyp,input,muXp(lags),SigX(lags,lags),post);
    m(k) = mtmp; 
    
    if s2tmp < 0
        warning(strcat([fun_name,': negative variance, step:',num2str(k)]))
    end

    s2(k) = s2tmp;
    
    [mutmp, sig2tmp] = gpx(hyp, inf, mean,cov,lik, input, target, muX(lags),post);
    
    mu(k) = mutmp;   
    
    if sig2tmp < 0
        warning(strcat([fun_name,': no propagation - negative variance, step:',...
            num2str(k)])); 
     
    end
    
    sig2(k) = sig2tmp; 
end 

% transform into column vectors 
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

end


function [M, S, V] = predictExactSEard(hyp,input,m,s,post)

[n, D] = size(input);    % number of examples and dimension of inputs
X = unwrap(hyp);                              % short hand for hyperparameters
k = zeros(n,1); M = 0; V = zeros(D,1); S = 0;
inp = bsxfun(@minus,input,m);                     % centralize inputs

invQ  =post.invQ;
alphaQ=post.alphaQ;
L	  =post.L;

% 2) compute predicted mean and inv(s) times input-output covariance
  
  iL = diag(exp(-X(1:D))); % inverse length-scales, inverse Lambda, iL
  in = inp*iL;
  B = iL*s*iL+eye(D); 
  
  t = in/B;
  l = exp(-sum(in.*t,2)/2);
  lb = l.*alphaQ;
  tiL = t*iL;
  c = exp(2*X(D+1))/sqrt(det(B));
  
  M = sum(lb)*c;                                           % predicted mean
  V = tiL'*lb*c;                    % inv(s) * input-output covariance
  k = 2*X(D+1)-sum(in.*in,2)/2;

% 3) compute predictive covariance, non-central moments
               
  ii = bsxfun(@rdivide,inp,exp(2*X(1:D)'));
    R = s*diag(exp(-2*X(1:D))+exp(-2*X(1:D)))+eye(D); 
    t = 1/sqrt(det(R));
    L = exp(bsxfun(@plus,k,k')+maha(ii,-ii,R\s/2));
    S = t*(alphaQ'*L*alphaQ - sum(sum(invQ.*L)));
  
  S = S + exp(2*X(D+1)) + exp(2*X(D+2));

% 4) centralize moments
S = S - M*M';                                              
end

function K = maha(a, b, Q)                         
  aQ = a*Q;
  K = bsxfun(@plus,sum(aQ.*a,2),sum(b*Q.*b,2)')-2*aQ*b';
end
