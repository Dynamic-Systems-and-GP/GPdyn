function [out1, out2, out3, out4] = gpSD00(hyp, inf, mean, cov, lik, input, target, ...
							targetvariance, derivinput, derivtarget, derivvariance, xs)
% gpSD00 computes the predictive mean and variance at test input for
% the GP model with incorporated local models, i.e., derivative observations, 
% with the covariance function covSEard and gaussian likelihood likGauss.
% It is used for training and prediction. 
%
%
%% Syntax
%  training:   [nlZ dnlZ] = gpSD00(hyp, inf, mean, cov, lik, input, target, 
%								 derivinput, derivtarget, derivvariance)
%
%  prediction: [mu S2 muderiv S2deriv]  = gpSD00(hyp,inf,mean,cov,lik,
%		           input, target, derivinput, derivtarget, derivvariance, xs)
%
%% Description
% This function computes the predictive mean and variance at test input for
% the LMGP model with the covariance function as sum of covSEard and
% covNoise. The form of the covariance function is
%
% C(x^p,x^q) = sf^2 * exp[-(x^p - x^q)'*inv(P)*(x^p - x^q)/2]
%            + sn^2 * delta_{p,q}
%
% where the first term is the squared negative exponential and the second term
% with the kronecker delta is the noise contribution. The P matrix is diagonal
% with "Automatic Relevance Determination" (ARD) or "input length scale"
% parameters ell_1^2,...,ell_D^2; The hyperparameter sf is the "signal std dev" and
% sn is the "noise std dev". All hyperparameters are collected in the vector X
% as follows:
%
% hyp.cov = [ log(ell_1)
%        	 ...
%       	 log(ell_D)
%       	 log(sf)]
%
% hyp.lik  =[ log(sn) ]  
% hyp.mean =[ ]
%
% Note: the reason why the log of the parameters are used in X is that this
% often leads to a better conditioned (and unconstrained) optimization problem
% than using the raw hyperparameters themselves.
%
% This function can conveniently be used with the "minimize" function to train
% a Gaussian process:
%
% [nlZ, dnlZ, i] = minimize(hyp, 'gpSD00', length, input, target)
%
%      
% Input: 
% * hyp            ... a struct of hyperparameters
%   inf      	   ... the inference method 	  --> this is never used here
%   cov      	   ... prior covariance function  --> this is never used here
%   mean    	   ... prior mean function        --> this is never used here
%   lik      	   ... likelihood function        --> this is never used here
% * input          ... a n by D matrix of training inputs
% * target         ... a (column) vector (of size n) of targets
% * targetvariance ... a (column) vector (of size n) of variances of
%                      target's unknown variances are indicated by NaN - these are then
%                      replaced by the noise term in the covariance function
% * derivinput     ... an n by D matrix of training inputs at which we have
%                      derivative information (not necessarily the same as 'input').
% * derivtarget    ... n by D matrix of partial derivatives at 'derivinput',
%                      w.r.t. each input
% * derivvariance  ... an n by D^2 matrix, where each row is the elements of
%                      the covariance matrix associated with the appropriate derivtarget
% * xs             ... a nn by D matrix of test inputs
%
% Output: 
% * nlZ            ... the returned value of negative log likelihood
% * dnlZ           ... a (column) vector (of size D+2) of partial derivatives
%                      of negative log likelihood w.r.t. the hyperparameters
% * mu             ... a (column) vector (of size nn) of predicted means
% * S2             ... a (column) vector (of size nn) of predicted variances (including noise variance)
% * muderiv        ... a vector of derivatives of predicted mean
% * S2deriv        ... derivatives of predicted variances
%
% where D is the dimension of the input. 



%
% Examples
% demo_example_lmgp_simulation.m
%
% See Also
% simullmgp00naive.m, simullmgp00mcmc.m, minimize.m
%
%% Signature
%  Copyright (C) by Carl Edward Rasmussen 1999 - 2003 (2003-07-24).
%    * modified 2003-07-25 by Roderick Murray-Smith. Derivative adaptation. Can now cope
%      with a number of derivative observations independent of the number of
%      function observations. It has seperate noise level for the derivative
%      observations.
%    * modified 2013 by Jus Kocijan. Comments.
%    * modified 2015 by Martin Stepancic. The input arguments match the gp.m script from gpml >= 3.1 toolbox

[n, D]      = size(input);              % number of examples and dimension of input space
[nD, D] = size(derivinput);             % number of derivative examples and dimension of input space
X=[hyp.cov(:);hyp.lik(:)];

if isempty(inf),  inf = @infExact; else                        % set default inf
  if iscell(inf), inf = inf{1}; end                      % cell input is allowed
  if ischar(inf), inf = str2func(inf); end        % convert into function handle
end
if isempty(mean),mean = {@meanZero}; end 			    % set default mean
if ischar(mean)  || isa(mean,  'function_handle'), mean  = {mean};  end  % make cell
mean1 = mean{1}; if isa(mean1, 'function_handle'), mean1 = func2str(mean1); end
if isempty(cov),cov = {@covSEard}; end 			    % set default covariance
if ischar(cov)  || isa(cov,  'function_handle'), cov  = {cov};  end  % make cell
cov1 = cov{1}; if isa(cov1, 'function_handle'), cov1 = func2str(cov1); end
if strcmp(cov1,'covFITC'); inf = @infFITC; end       % only one possible inf alg
if isempty(lik),  lik = @likGauss; else                        % set default lik
  if iscell(lik), lik = lik{1}; end                      % cell input is allowed
  if ischar(lik), lik = str2func(lik); end        % convert into function handle
end

if strcmpi(cov1,'covSEard')==0 error('Unsupported covariance function. Please use covSEard.m'); end
if strcmpi(func2str(lik),'likGauss')==0 error('Unsupported likelihood function. Please use likGauss.m'); end
if strcmpi(mean1,'meanZero')==0 error('Unsupported mean function. Please use meanZero.m'); end
if (size(X,1)~=D+2) error('Number of hyperparameters disagree with input dataset dimensions'); end

% create the full input matrix, stacking the function observations and
% derivative observations (one repeat for each partial derivative) together
N = n+D*nD;
fullinput  = [input; repmat(derivinput,D,1)];
input      = input;%./repmat(exp(X(1:D))',n,1);
fullinput  = fullinput;% ./ repmat(exp(X(1:D))',N,1);
derivinput = derivinput;%./ repmat(exp(X(1:D))',nD,1);

fulltarget = [target; derivtarget(:)];

% first, we write out the covariance matrix Q. This includes covariance of
% all inputs, and derivative point inputs (N = n+D*nD). The first n x n
% matrix corresponds to the classical GP cov. matrix
Z = zeros(N,N);
for d = 1:D
  Z = Z + (repmat(fullinput(:,d),1,N)-repmat(fullinput(:,d)',N,1)).^2*exp(X(d));
end
Z = exp(2*X(D+1))*exp(-0.5*Z);

Q=Z;

% now fill in the blocks corresponding to the covariance between
% derivatives and function observations, and among the various partial
% derivatives
for d=1:D
    % distance among derivative inputs
    Z1 = (repmat(derivinput(:,d),1,nD)-repmat(derivinput(:,d)',nD,1))*exp(X(d));
    % distance between derivative inputs and function inputs
    Z2 = (repmat(input(:,d),1,nD)-repmat(derivinput(:,d)',n,1))*exp(X(d));
    
    % calculate cross-terms (function-deriv) in cov. matrix)
    Q(1:n, n+(d-1)*nD+1:d*nD+n) = Q(1:n, n+(d-1)*nD+1:d*nD+n).*Z2;
    Q(n+(d-1)*nD+1:d*nD+n,1:n)  = Q(n+(d-1)*nD+1:d*nD+n,1:n).*Z2';
    
    for j=1:D
        % calculate covariance among different derivative terms
        Q(n+(j-1)*nD+1:j*nD+n,n+(d-1)*nD+1:d*nD+n)  =  Q(n+(j-1)*nD+1:j*nD+n,n+(d-1)*nD+1:d*nD+n).*Z1;
        Q(n+(d-1)*nD+1:d*nD+n, n+(j-1)*nD+1:j*nD+n) =  Q(n+(d-1)*nD+1:d*nD+n, n+(j-1)*nD+1:j*nD+n).*Z1';
    end
    % calculate covariance among same derivative terms
    Q(n+(d-1)*nD+1:d*nD+n,n+(d-1)*nD+1:d*nD+n) = Q(n+(d-1)*nD+1:d*nD+n,n+(d-1)*nD+1:d*nD+n) ...
        +Z(n+(d-1)*nD+1:d*nD+n,n+(d-1)*nD+1:d*nD+n)*exp(X(d));
end

% reformat the covariance matrices of the derivative estimates from columns
% to matrices, spread appropriately using kronecker products.
derivcovar = zeros(nD*D,nD*D);
for i=1:nD
    pos = zeros(nD,nD); 
    pos(i,i) = 1;
    derivcovar = derivcovar + kron(reshape(derivvariance(i,:),D,D),pos);
end

unknownvarind = find(isnan(targetvariance));
knownvarind = find(isfinite(targetvariance));
noisediag = zeros(N,N);
noisediag(unknownvarind,unknownvarind) = diag(repmat(exp(2*X(D+2)),length(unknownvarind),1));  % first points have no known variance
noisediag(knownvarind,knownvarind) = diag(targetvariance(knownvarind));                   % then points with known variance
noisediag(n+1:end,n+1:end) = derivcovar; %exp(2*derivvariance(:))]);     % then derivative points with known variance
noisediag = noisediag + 1e-5*eye(N,N); %jitter

if nargin == 11   % if no test cases, we compute the negative log likelihood ...

  W = inv(Q+noisediag);               % W is inv (Q plus noise term)
  invQt = W*fulltarget;                               % don't compute determinant..
  logdetQ = 2*sum(log(diag(chol(Q+noisediag))));        % ..directly
  out1 = 0.5*logdetQ + 0.5*fulltarget'*invQt + 0.5*N*(D+1)*log(2*pi);

  % ... and its partial derivatives

  out2 = zeros(D+2,1);                  % set the size of the derivative vector
  W = W-invQt*invQt';
  QW = W.*Q;
  for d = 1:D   % derivative w.r.t. d-th input
      dist=-0.5*(repmat(fullinput(:,d),1,N)-repmat(fullinput(:,d)',N,1)).^2*exp(X(d));
      % deriv of cov. of func. observations
      V= Q.*dist;
      % now add in cross terms
      V(1:n,n+(d-1)*nD+1:d*nD+n) = V(1:n,n+(d-1)*nD+1:d*nD+n)+Q(1:n,n+(d-1)*nD+1:d*nD+n);
      V(n+(d-1)*nD+1:d*nD+n,1:n) = V(n+(d-1)*nD+1:d*nD+n,1:n)+Q(n+(d-1)*nD+1:d*nD+n,1:n);
     
      % derivative covariances.
      for j=1:D
          V(n+(j-1)*nD+1:j*nD+n,n+(d-1)*nD+1:d*nD+n) = V(n+(j-1)*nD+1:j*nD+n,n+(d-1)*nD+1:d*nD+n)+Q(n+(j-1)*nD+1:j*nD+n,n+(d-1)*nD+1:d*nD+n);
          V(n+(d-1)*nD+1:d*nD+n,n+(j-1)*nD+1:j*nD+n) = V(n+(d-1)*nD+1:d*nD+n,n+(j-1)*nD+1:j*nD+n)+Q(n+(d-1)*nD+1:d*nD+n,n+(j-1)*nD+1:j*nD+n);
      end
      V(n+(d-1)*nD+1:d*nD+n,n+(d-1)*nD+1:d*nD+n) = V(n+(d-1)*nD+1:d*nD+n,n+(d-1)*nD+1:d*nD+n) ...
                                    -exp(X(d)).*Z(n+(d-1)*nD+1:d*nD+n,n+(d-1)*nD+1:d*nD+n);
      out2(d) = (-invQt'*V*invQt+sum(sum(inv(Q+noisediag).*V)))/2;

  end
  out2(D+1) = sum(sum(QW));
  out2(D+2) = trace(W(unknownvarind,unknownvarind))*exp(2*X(D+2));

else                    % ... otherwise compute (marginal) test predictions ...
  [nn, D] = size(xs);     % number of test cases and dimension of input space
  xs = xs;% ./ repmat(exp(X(1:D))',nn,1);

  a = zeros(N, nn);    % compute the covariance between training and test cases
  for d = 1:D
    a = a + (repmat(fullinput(:,d),1,nn)-repmat(xs(:,d)',N,1)).^2*exp(X(d));
  end
  a = exp(2*X(D+1))*exp(-0.5*a);
  Z = a;
  for d=1:D
     ZZ = (repmat(derivinput(:,d),1,nn)-repmat(xs(:,d)',nD,1))*exp(X(d));
     a(n+(d-1)*nD+1:d*nD+n,1:nn) =  -a(n+(d-1)*nD+1:d*nD+n,1:nn) .* ZZ;     
  end

  % ... write out the desired terms

  if nargout == 1
      out1 = a'*((Q+noisediag)\fulltarget);              % predicted means
  else
      invQ = inv(Q+noisediag);
      out1 = a'*(invQ*fulltarget);                       % predicted means
      out2 = exp(2*X(D+2)) + exp(2*X(D+1)) - sum(a.*(invQ*a),1)';								  % predicted variance (including noise variance)
      for d = 1:D
          c = a .* (repmat(fullinput(:,d),1,nn)-repmat(xs(:,d)',N,1));
          c(n+(d-1)*nD+1:d*nD+n,1:nn) = c(n+(d-1)*nD+1:d*nD+n,1:nn) + Z(n+(d-1)*nD+1:d*nD+n,1:nn);
          out3(1:nn,d) = exp(X(d))*c'*invQ*fulltarget;                    % derivative of mean
          out4(1:nn,d) = exp(X(d))*(exp(2*X(D+1))-exp(X(d))*sum(c.*(invQ*c),1)');
      end
  end
end
