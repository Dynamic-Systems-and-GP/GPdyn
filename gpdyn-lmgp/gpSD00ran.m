%% gpSD00ran
function [mu, S2, invQ] = gpSD00ran(X, input, target, targetvariance, derivinput, derivtarget, derivvariance, test, SigmaX)

%% Syntax
%  function [mu, S2, invQ] = gpSD00ran(X, input, target, targetvariance, derivinput, derivtarget, derivvariance, test, SigmaX)

%% Description
% This function computes the predictive mean and variance at test input for
% the LMGP model with the covariance function as sum of covSEard and
% covNoise.  The prediction and propagation of variance is done by this
% routine. When test data are given, then (marginal) Gaussian predictions
% are computed, whose mean and (noise free) variance are returned.
% 
% Inputs: 
% X      is a (column) vector (of size D+2) of hyperparameters
% input  is a n by D matrix of training inputs
% target is a (column) vector (of size n) of targets
% targetvariance is a (column) vector (of size n) of variances of
% targets. unknown variances are indicated by NaN - these are then
%   replaced by the noise term in the covariance function
% derivinput is an n by D matrix of training inputs at which we have
%   derivative information (not necessarily the same as 'input').
% derivtarget is an n by D matrix of partial derivatives at 'derivinput',
%   w.r.t. each input
% derivvariance is an n by D^2 matrix, where each row is the elements of
%   the covariance matrix associated with the appropriate derivtarget
% test   is a nn by D matrix of test inputs
% SigmaX   covariance of the test input 
% Outputs: 
% mu     is a (column) vector (of size nn) of prediced means
% S2     is a (column) vector (of size nn) of predicted variances
% invQ   is the inverse or covariance matrix.
%      
% (C) Copyright 1999 - 2003, Carl Edward Rasmussen (2003-07-24).
% derivative adaptation, Roderick Murray-Smith (2003-07-25). Can now cope
% with a number of derivative observations independent of the number of
% function observations. It has seperate noise level for the derivative
% observations.
% Variance propagation upgrade, Jus Kocijan, 2003

%% Examples
% demo_example_lmgp_simulation.m

%% See Also
% SIMULLMGP00NAIVE SIMULLMGP00MCMC, MINIMIZE, GPSD00


[n, D]      = size(input);              % number of examples and dimension of input space
[nD, D] = size(derivinput);             % number of derivative examples and dimension of input space


% create the full input matrix, stacking the function observations and
% derivative observations (one repeat for each partial derivative) together
N = n+D*nD;
fullinput  = [input; repmat(derivinput,D,1)];
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

% computation of (marginal) test predictions with propagation

invQ = inv(Q+noisediag); % total inverse matrix

[nn, D] = size(test);     % number of test cases and dimension of input space

 W=inv(diag(exp(X(1:D))));
 a = zeros(N, nn);    % compute the covariance between training and test cases
 V=inv(W+SigmaX);
  for d=1:D
          for e = 1:D
          a = a + (repmat(fullinput(:,d),1,nn)-repmat(test(:,e)',N,1)).^2*V(d,e); 
          end
        end


  a = exp(2*X(D+1))*det(eye(D)+inv(W)*SigmaX)^(-0.5)*exp(-0.5*a);
      
  % ... calculation of mean value for propagation    
     
     c=repmat(derivinput,D,1)*(eye(D)-inv(W+SigmaX)*W)+repmat(test,D*nD,1)*inv(inv(W)*SigmaX+eye(D));
     ad=[];
     for d=1:D
         ad=[ad;
            -exp(X(d)).*a(n+(d-1)*nD+1:d*nD+n,1:nn).*(fullinput(n+(d-1)*nD+1:d*nD+n,d)-c((d-1)*nD+1:d*nD,d))];
        end
     a=[a(1:n,1);
        ad];
    
     out1=a'*(invQ*fulltarget);
     
  % ... calculation of variance for propagation    

     beta=invQ*fulltarget;

  % compute the covariance Cmod2
     a2 = zeros(N,N);   
     for d = 1:D
          a2 = a2 + (repmat(fullinput(:,d),1,N)-repmat(fullinput(:,d)',N,1)).^2*inv(2*W(d,d));  
     end

      a2 = exp(2*X(D+1))*2^(-D*0.5)*exp(-0.5*a2);
  
    % compute the covariance Cmod3
     a3 = zeros(N, N);
     V=inv(W/2+SigmaX);
     for d=1:D
          for e = 1:D
          a3 = a3 + (repmat(test(:,d)',N,N)-(repmat(fullinput(:,e),1,N)+repmat(fullinput(:,e)',N,1))/2).^2*V(d,e); 
          end
        end
 
     a3 = exp(2*X(D+1))*det(0.5*eye(D)+inv(W)*SigmaX)^(-0.5)*exp(-0.5*a3);

% compute the CX matrix 
     
     C=W/2-W/2*inv(W/2+SigmaX)*W/2;
     help1=eye(D)-inv(W/2+SigmaX)*W/2;
     help2=test*inv(inv(W/2)*SigmaX+eye(D));
   
     cij=ones(N,N,D);
     CX=ones(N,N);  
     for d=1:D
         for i=1:N
            for j=1:N
                cij(i,j,d)=0.5*(fullinput(i,:)+fullinput(j,:))*help1(:,d)+help2(:,d);
            end
        end

     
        for i=1:n
            for j=n+(d-1)*nD+1:d*nD+n
                CX(i,j)=-inv(W(d,d))*(fullinput(j,d)-cij(i,j,d));
            end
        end
        for i=n+(d-1)*nD+1:d*nD+n
            for j=1:n
                CX(i,j)=-inv(W(d,d))*(fullinput(i,d)-cij(i,j,d));
            end
        end
    end
   
    for di=1:D
        for dj=1:D
            for i=n+(di-1)*nD+1:di*nD+n
                for j=n+(dj-1)*nD+1:dj*nD+n
                CX(i,j)=inv(W(di,di))*inv(W(dj,dj))*(fullinput(i,di)*fullinput(j,dj)-(fullinput(j,dj)*cij(i,j,di)+fullinput(i,di)*cij(i,j,dj))+C(di,dj)+cij(i,j,di)*cij(i,j,dj));
                end    
            end
        end
    end   

    s2d=sum(sum((invQ-beta*beta').*CX.*a2.*a3));
    out2=exp(2*X(D+1)) - s2d - out1^2;
    
    mu=out1;
    S2=out2;
