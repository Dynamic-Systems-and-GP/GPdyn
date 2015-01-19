%% loss 
function [ae, se, lpd, mrse] = loss2(tt, mu, s2)

%% Syntax
% function [ae, se, lpd, mrse] = loss2(tt, mu, s2)


%% Description
% Function computes several frequently used performance values: 
% - mean absolute error AE 
% - mean squared error SE 
% - minus log-predicted density error LPD 
% -  mean relative squared error MRSE 
% Notes: 
% (1) It can be used with one-step-ahead prediction or simulation results. 
% (2) Can be expanded with other performance values. 
% Inputs: 
% tt .. target output values (=system output) 
% y .. predicted output values (=model output) 
% s2 .. predicted output variance 
% Outputs: 
% ae .. mean absolute error 
% se .. mean squared error 
% lpd .. minus log-predicted density error 
% mrse .. mean relative square error 


if size(tt) ~= size(mu)
  mu = mu';  
end


ae = mean(abs(tt-mu));
se = mean((tt-mu).^2);

% mean relative squared error MRSE 
mrse = sqrt(sum((tt-mu).^2)/sum(tt.^2));

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

lpd = 0.5*mean(log(2*pi) + log(s2) + (tt-mu).^2./s2);
printout(ae,se);
printoutlpd(lpd);
printoutmrse(mrse);
return; 


function []=printout(ae,se);
disp(strcat('AE = ',num2str(ae)));
disp(strcat('SE = ',num2str(se)));

function []=printoutlpd(lpd);
disp(strcat('LD = ',num2str(lpd)));

function []=printoutmrse(mrse);
disp(strcat('MRSE = ',num2str(mrse)));



