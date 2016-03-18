function K = covMaternard(d, hyp, x, z, i)

% Matern covariance function with nu = d/2 and Automatic Relevance Detemination
% (ARD) distance measure. For % d=1 the function is also known as the exponential 
% covariance function or the Ornstein-Uhlenbeck covariance. The covariance 
% function is:
%
%   k(x^p,x^q) = s2f * f( sqrt(d)*r ) * exp(-sqrt(d)*r)
%
% with f(t)=1 for d=1, f(t)=1+t for d=3 and f(t)=1+t.*(1+t/3) for d=5.
% Here r is the distance sqrt((x^p-x^q)'*inv(P)*(x^p-x^q)), P is ell times
% the unit matrix and sf2 is the signal variance. The hyperparameters are:
%
% hyp = [ log(ell_1)
%         log(ell_2)
%          .
%         log(ell_D)
%         log(sqrt(sf2)) ]
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-09-10.
%                  Ali Abusnina, 2014
%
% See also COVFUNCTIONS.M.

if nargin<3, K = '(D+1)'; return; end              % report number of parameters
if nargin<4, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

[n,D] = size(x);
ell = exp(hyp(1:D));                               % characteristic length scale
sf2 = exp(2*hyp(D+1));  

switch d
  case 1, f = @(t) 1;               df = @(t) 1;
  case 3, f = @(t) 1 + t;           df = @(t) t;
  case 5, f = @(t) 1 + t.*(1+t/3);  df = @(t) t.*(1+t)/3;
end
          m = @(t,f) f(t).*exp(-t); dm = @(t,f) df(t).*t.*exp(-t);

% precompute distances
if dg                                                               % vector kxx
  K = zeros(size(x,1),1);
else
  if xeqz                                                 % symmetric matrix Kxx
      K = sqrt( sq_dist(sqrt(d)* (bsxfun(@rdivide, x', ell))    ));

  else                                                   % cross covariances Kxz
    
    K = sqrt( sq_dist(sqrt(d)* (bsxfun(@rdivide, x', ell)) , sqrt(d) * (bsxfun(@rdivide, z', ell))) );
    
  end
end

if nargin<5                                                        % covariances
  K = sf2*m(K,f);
else                                                               % derivatives
  if i<=D 
    K = sf2*dm(K,f);
  elseif i<=D+1
    K = 2*sf2*m(K,f);
  else
    error('Unknown hyperparameter')
  end
end