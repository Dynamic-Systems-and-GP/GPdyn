function [out1, out2, out3, out4] = gpSD00(hyp, input, target, targetvariance, derivinput, derivtarget, derivvariance, test)
% gpSD00 computes the predictive mean and variance at test input for
% the GP model with incorporated local models, i.e., derivative observations, 
% with the covariance function as sum of covSEard and
% covNoise. It is used for training and prediction. 
%
%
%% Syntax
%  training:   [fX dfX] = gpSD00(X, input, target, derivinput, derivtarget, derivvariance)
%  prediction: [mu S2 muderiv S2deriv]  = gpSD00(X, input, target, derivinput,
% derivtarget, derivvariance, test)
%
%% Description
% This function computes the predictive mean and variance at test input for
% the LMGP model with the covariance function as sum of covSEard and
% covNoise. The form of the covariance function is
%
% C(x^p,x^q) = v^2 * exp[-(x^p - x^q)'*inv(P)*(x^p - x^q)/2]
%            + u^2 * delta_{p,q}
%
% where the first term is the squared negative exponential and the second term
% with the kronecker delta is the noise contribution. The P matrix is diagonal
% with "Automatic Relevance Determination" (ARD) or "input length scale"
% parameters l_1^2,...,l_D^2; The hyperparameter v is the "signal std dev" and
% u is the "noise std dev". All hyperparameters are collected in the struct X.
%
% Note: the reason why the log of the parameters are used in X is that this
% often leads to a better conditioned (and unconstrained) optimization problem
% than using the raw hyperparameters themselves.
%
% This function can conveniently be used with the "minimize" function to train
% a Gaussian process.
%
%
%   training: [nlZ dnlZ          ] = gp(hyp, inf, mean, cov, lik, x, y);
% prediction: [ymu ys2 fmu fs2   ] = gp(hyp, inf, mean, cov, lik, x, y, xs);
%         or: [ymu ys2 fmu fs2 lp] = gp(hyp, inf, mean, cov, lik, x, y, xs, ys);
%
% where:
%
%   hyp      column vector of hyperparameters
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
% * test           ... a nn by D matrix of test inputs
%
% Output: 
% * fX             ... the returned value of minus log likelihood
% * dfX            ... a (column) vector (of size D+2) of partial derivatives
%                      of minus the log likelihood wrt each of the hyperparameters
% * mu             ... a (column) vector (of size nn) of predicted means
% * S2             ... a (column) vector (of size nn) of predicted variances
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
%    * modified 2015 by Martin Stepancic

[n, D]      = size(input);              % number of examples and dimension of input space
[nD, D] = size(derivinput);             % number of derivative examples and dimension of input space


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
  Z = Z + (repmat(fullinput(:,d),1,N)-repmat(fullinput(:,d)',N,1)).^2*exp(hyp(d));
end
Z = exp(2*hyp(D+1))*exp(-0.5*Z);

Q=Z;

% now fill in the blocks corresponding to the covariance between
% derivatives and function observations, and among the various partial
% derivatives
for d=1:D
    % distance among derivative inputs
    Z1 = (repmat(derivinput(:,d),1,nD)-repmat(derivinput(:,d)',nD,1))*exp(hyp(d));
    % distance between derivative inputs and function inputs
    Z2 = (repmat(input(:,d),1,nD)-repmat(derivinput(:,d)',n,1))*exp(hyp(d));
    
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
        +Z(n+(d-1)*nD+1:d*nD+n,n+(d-1)*nD+1:d*nD+n)*exp(hyp(d));
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
noisediag(unknownvarind,unknownvarind) = diag(repmat(exp(2*hyp(D+2)),length(unknownvarind),1));  % first points have no known variance
noisediag(knownvarind,knownvarind) = diag(targetvariance(knownvarind));                   % then points with known variance
noisediag(n+1:end,n+1:end) = derivcovar; %exp(2*derivvariance(:))]);     % then derivative points with known variance
noisediag = noisediag + 1e-5*eye(N,N); %jitter

if nargin == 7   % if no test cases, we compute the negative log likelihood ...

  W = inv(Q+noisediag);               % W is inv (Q plus noise term)
  invQt = W*fulltarget;                               % don't compute determinant..
  logdetQ = 2*sum(log(diag(chol(Q+noisediag))));        % ..directly
  out1 = 0.5*logdetQ + 0.5*fulltarget'*invQt + 0.5*N*(D+1)*log(2*pi);

  % ... and its partial derivatives

  out2 = zeros(D+2,1);                  % set the size of the derivative vector
  W = W-invQt*invQt';
  QW = W.*Q;
  for d = 1:D   % derivative w.r.t. d-th input
      dist=-0.5*(repmat(fullinput(:,d),1,N)-repmat(fullinput(:,d)',N,1)).^2*exp(hyp(d));
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
                                    -exp(hyp(d)).*Z(n+(d-1)*nD+1:d*nD+n,n+(d-1)*nD+1:d*nD+n);
      out2(d) = (-invQt'*V*invQt+sum(sum(inv(Q+noisediag).*V)))/2;

  end
  out2(D+1) = sum(sum(QW));
  out2(D+2) = trace(W(unknownvarind,unknownvarind))*exp(2*hyp(D+2));

else                    % ... otherwise compute (marginal) test predictions ...

  [nn, D] = size(test);     % number of test cases and dimension of input space
  test = test;% ./ repmat(exp(X(1:D))',nn,1);

  a = zeros(N, nn);    % compute the covariance between training and test cases
  for d = 1:D
    a = a + (repmat(fullinput(:,d),1,nn)-repmat(test(:,d)',N,1)).^2*exp(hyp(d));
  end
  a = exp(2*hyp(D+1))*exp(-0.5*a);
  Z = a;
  for d=1:D
     ZZ = (repmat(derivinput(:,d),1,nn)-repmat(test(:,d)',nD,1))*exp(hyp(d));
     a(n+(d-1)*nD+1:d*nD+n,1:nn) =  -a(n+(d-1)*nD+1:d*nD+n,1:nn) .* ZZ;     
  end

  % ... write out the desired terms

  if nargout == 1
      out1 = a'*((Q+noisediag)\fulltarget);              % predicted means
  else
      invQ = inv(Q+noisediag);
      out1 = a'*(invQ*fulltarget);                       % predicted means
      out2 = exp(2*hyp(D+1)) - sum(a.*(invQ*a),1)'; % predicted noise-free variance
      for d = 1:D
          c = a .* (repmat(fullinput(:,d),1,nn)-repmat(test(:,d)',N,1));
          c(n+(d-1)*nD+1:d*nD+n,1:nn) = c(n+(d-1)*nD+1:d*nD+n,1:nn) + Z(n+(d-1)*nD+1:d*nD+n,1:nn);
          out3(1:nn,d) = exp(hyp(d))*c'*invQ*fulltarget;                    % derivative of mean
          out4(1:nn,d) = exp(hyp(d))*(exp(2*hyp(D+1))-exp(hyp(d))*sum(c.*(invQ*c),1)');
      end
  end
end
