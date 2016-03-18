function [ae, se, lpd, mrse, smse, msll] = loss(y, yp_m, yp_se2)
% Function computes several frequently used performance measures
%
%% Syntax
% function [ae, se, lpd, mrse, smse, msll] = loss(y, yp_m, yp_se2)
%
%
%% Description
% Function computes several frequently used performance measures. It can be
% used with one-step-ahead prediction or simulation results. 
%
% Input: 
% * y      ... target output values (=system output) 
% * yp_m   ... predicted output mean values (=model output) 
% * yp_se2 ... predicted output variances   (=sigma^2)
% 
% Output: 
% * ae     ... the mean absolute error 
% * se     ... the mean squared error 
% * lpd    ... the log-predicted loss
% * mrse   ... the mean relative square error 
% * smse   ... the standardized mean squared error
% * msll   ... the mean standardized log loss
% Examples:
% demo_example_gp_simulation.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Verification of input arguments
if nargin < 2  error('not enough input arguments'); end
y=y(:);
yp_m = yp_m(:);
if size(y) ~= size(yp_m)  error('wrong sizes'); end
if nargin == 3
  yp_se2 = yp_se2(:);
  if size(y) ~= size(yp_se2) error('wrong sizes');  end
else
  yp_se2=[];
end
%%%

y_m=mean(y);    % mean of output
y_se2=var(y);   % variance of output

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AE: (mean absolute error )
ae = mean(abs(y-yp_m));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SE: (mean squared error)
se = mean((y-yp_m).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MRSE: mean relative squared error
mrse = sqrt(sum((y-yp_m).^2)/sum(y.^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SMSE: standardized mean squared error
smse = mean((yp_m-y).^2)/(y_se2);

if ~isempty(yp_se2)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % LPD: log-predicted density loss
 lpd =      0.5*mean(log(2*pi) + log(yp_se2) + (y-yp_m).^2./yp_se2);  
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % MSLL: mean standardized log loss
 msll = lpd-0.5*mean(log(2*pi) + log(y_se2 )  + (y-y_m ).^2./y_se2 );
else
 lpd =  'undefined, no variance provided';
 msll = 'undefined, no variance provided';
end

printout('AE  ',ae);
printout('SE  ',se);
printout('MRSE',mrse);
printout('SMSE',smse);
printout('LPD ',lpd);
printout('MSLL',msll);

return; 

function []=printout(errtype,value)
disp([errtype ' = ',num2str(value)]);



