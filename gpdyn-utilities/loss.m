function [ae, se, lpd, mrse] = loss(tt, y, s2)
% Function computes several frequently used perforamnce values
%
%% Syntax
% function [ae, se, lpd, mrse] = loss2(tt, y, s2)
%
%
%% Description
% Function computes several frequently used perforamnce values. It can be
% used with one-step-ahead prediction or simulation results. 
% * mean absolute error AE 
% * mean squared error SE 
% * minus log-predicted density error LPD 
% * mean relative squared error MRSE 
%
% Input: 
% * tt   ... target output values (=system output) 
% * y    ... predicted output values (=model output) 
% * s2   ... predicted output variance 
% 
% Output: 
% * ae   ... mean absolute error 
% * se   ... mean squared error 
% * lpd  ... minus log-predicted density error 
% * mrse ... mean relative square error 
%
% Examples:
% demo_example_gp_simulation.m


if size(tt) ~= size(y)
  y = y';  
end


ae = mean(abs(tt-y));
se = mean((tt-y).^2);

% mean relative squared error MRSE 
mrse = sqrt(sum((tt-y).^2)/sum(tt.^2));

if(nargin==2)
    lpd = []; 
    disp('loss: input s2 is undefined');
    printout(ae,se);
    printoutmrse(mrse);
    return; 
elseif (isempty(s2))
    lpd = [];
    printout(ae,se);
    printoutmrse(mrse);
    return; 
end

if size(tt) ~= size(s2)  
  s2 = s2';
end

lpd = 0.5*mean(log(2*pi) + log(s2) + (tt-y).^2./s2);
printout(ae,se);
printoutlpd(lpd);
printoutmrse(mrse);
return; 


function []=printout(ae,se);
disp(strcat('AE = ',num2str(ae)));
disp(strcat('SE = ',num2str(se)));

function []=printoutlpd(lpd);
disp(strcat('LPD = ',num2str(lpd)));

function []=printoutmrse(mrse);
disp(strcat('MRSE = ',num2str(mrse)));



