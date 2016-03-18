function K = covNNard(hyp, x, z, i)

% Neural network covariance function with a single parameter for the distance
% measure. The covariance function is parameterized as:
%
% k(x^p,x^q) = sf2 * asin(x^p'*P*x^q / sqrt[(1+x^p'*P*x^p)*(1+x^q'*P*x^q)])
%
% where the x^p and x^q vectors on the right hand side have an added extra bias
% entry with unit value. P is ell^-2 times the unit matrix and sf2 controls the
% signal variance. The hyperparameters are:
%
% hyp = [ log(ell)
%         log(sqrt(sf2) ]
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2010-09-10.
%
% See also COVFUNCTIONS.M.

if nargin<2, K = 'D+2'; return; end                  % report number of parameters
if nargin<3, z = []; end                                   % make sure, z exists
xeqz = numel(z)==0; dg = strcmp(z,'diag') && numel(z)>0;        % determine mode

[n,D] = size(x);
ell2 = exp(2*hyp(1:D));                               % characteristic length scale
ellb2 = exp(2*hyp(D+1));							  % characteristic length scale for bias element
sf2 = exp(2*hyp(D+2));                                         % signal variance

if xeqz
    deps=diag(ones(n,1)*eps);
else
    deps=0;
end
XSX = 1/ellb2 + x./repmat(ell2',[n,1])*x';
if dg                                                               % vector kxx
  K0 = diag(XSX)./(diag(XSX)+1);
else
  if xeqz                                                 % symmetric matrix Kxx
    K0 = XSX./(sqrt(1+diag(XSX))*sqrt(1+diag(XSX))');
  else                                                   % cross covariances Kxz
    nz = size(z,1);
    XSZ = 1/ellb2 + (x./repmat(ell2',[n,1]))*z';
    ZSZ = 1/ellb2 + (z./repmat(ell2',[nz,1]))*z';
    K0 = XSZ./(sqrt(1+diag(XSX))*sqrt(1+diag(ZSZ))');
  end
end

if nargin<4                                                        % covariances
  K = sf2*asin(K0-deps);
else                                                               % derivatives
  if i<=D                                                          % lengthscales
    sS = (x(:,i).^2)./(1+diag(XSX));
    if dg
        K =  -2*sf2/ell2(i)*( (-2*x(:,i).*x(:,i))./(1+diag(XSX)) - ...
             K0.*sS )./sqrt(1-K0.*K0);
    else
      if xeqz
        K =  -2*sf2/ell2(i)*( (x(:,i)*x(:,i)')./(sqrt(1+diag(XSX))*sqrt(1+diag(XSX))') - ...
             K0.*0.5.*(repmat(sS,[1,n]) + repmat(sS',[n,1])) )./sqrt(1-K0.*K0); % last division is asin derivative
      else  
        sZ = (z(:,i).^2)./(1+diag(ZSZ));
        K =  -2*sf2/ell2(i)*( (x(:,i)*z(:,i)')./(sqrt(1+diag(XSX))*sqrt(1+diag(ZSZ))') - ...
             K0.*0.5.*(repmat(sS,[1,nz]) + repmat(sZ',[n,1])) )./sqrt(1-K0.*K0); % last division is asin derivative
      end
      
    end   

  elseif i==D+1                                                        % lengthscale of bias
    sS = 1./(1+diag(XSX));
    if dg
        K =  -2*sf2/ellb2*(1./(1+diag(XSX)) - K0.*sS)./sqrt(1-K0.*K0);
    else
      if xeqz
        K =  -2*sf2/ellb2*( 1./(sqrt(1+diag(XSX))*sqrt(1+diag(XSX))') - ...
             K0.*0.5.*(repmat(sS,[1,n]) + repmat(sS',[n,1])) )./sqrt(1-K0.*K0); % last division is asin derivative
      else  
        sZ = 1./(1+diag(ZSZ));
        K =  -2*sf2/ellb2*( 1./(sqrt(1+diag(XSX))*sqrt(1+diag(ZSZ))') - ...
             K0.*0.5.*(repmat(sS,[1,nz]) + repmat(sZ',[n,1])) )./sqrt(1-K0.*K0); % last division is asin derivative
      end
    end   
  elseif i==D+2                                                        % magnitude
    K = 2*sf2*asin(K0-deps);
  else
    error('Unknown hyperparameter')
  end
end
