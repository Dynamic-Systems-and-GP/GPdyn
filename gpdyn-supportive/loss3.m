function [ae, se, lpd, mrse, smse, msll] = loss(tt, mu, s2)
% error measures
%
%% Syntax
% function [ae, se, lpd, mrse, smse, msll] = loss(tt, mu, s2)
%
%% Description
% Function computes several frequently used performance values: 
% - mean absolute error AE 
% - mean squared error SE 
% - minus log-predicted density error LPD 
% - mean relative squared error MRSE
% - standardized mean squared error MSE
% - mean standardized log loss MSLL
% 
% Notes: 
% (1) It can be used with one-step-ahead prediction or simulation results. 
% (2) Can be expanded with other performance values. 
% 
% Input:
% * tt ... target output values (=system output) 
% * y  ... predicted output values (=model output) 
% * s2 ... predicted output variance 
% 
% Output:
% * ae   .. mean absolute error 
% * se   .. mean squared error 
% * lpd  .. minus log-predicted density error 
% * mrse .. mean relative square error 
% * mse  .. standardized mean squared error
% * msll .. mean standardized log loss

if nargin < 2
  error('not enough input arguments');
else
  if size(tt) ~= size(mu)
    mu = mu';
    if size(tt) ~= size(mu)
      error('wrong sizes');
    end
  end
  if nargin == 3
    if size(tt) ~= size(s2)
      s2 = s2';
      if size(tt) ~= size(s2)
        error('wrong sizes');
      end
    end
  end
end

mtt = mean(tt);
vtt = var(tt);

% AE
ae = mean(abs(tt-mu));
%disp(strcat('AE   = ', num2str(ae)));
fprintf('AE   = %f\n', ae);

% SE
se = mean((tt-mu).^2);
%disp(strcat('SE   = ', num2str(se)));
fprintf('SE   = %f\n', se);

% MRSE
mrse = sqrt(sum((tt-mu).^2)/sum(tt.^2));
%disp(strcat('MRSE = ', num2str(mrse)));
fprintf('MRSE = %f\n', mrse);

% SMSE
smse = mean((mu-tt).^2)/(vtt);
%disp(strcat('SMSE  = ', num2str(smse)));
fprintf('SMSE  = %f\n', smse);

if nargin == 3
  % LPD
  lpd = 0.5*mean(log(2*pi) + log(s2) + (tt-mu).^2./s2);
  %disp(strcat('LD   = ', num2str(lpd)));
  fprintf('LD   = %f\n', lpd);
  
%   % MSLL
%   msll = 0.5*mean(((tt - mu)./sqrt(s2)).^2 + log(s2)) - 0.5*(vtt + mtt^2);
%   %disp(strcat('MSLL = ', num2str(msll)));
%   fprintf('MSLL = %f\n', msll);

  % MSLL
  msll = lpd-0.5*mean(log(2*pi) + log(vtt)+(tt-mtt).^2/vtt);
  %disp(strcat('MSLL = ', num2str(msll)));
  fprintf('MSLL = %f\n', msll);

end