function [MOIndex] = lipschitz(u, y, maxlag, model, fig)
% A method to determine the lag space, based on Lipschitz quotients
%
%% Syntax
% [MOIndex] = lipschitz(u, y, maxlag)
% 
%% Description
%  Given a set of corresponding inputs and outputs the function calculates
%  so called Lipschitz number for each combination of m and l where m
%  represents the number of delayed outputs, and l the number of delayed
%  inputs for the case of dynamic system: y(t) = f(y(t-1),...,y(t-l),
%  u(t-1),..., u(t-m)). 
%  To small m and l result in a large Lipschitz nuber, while to large lag
%  spaces do not have greater effect on the Lipschitz nuber. In order to 
%  determine the proper lag space we have to look for the knee point, where
%  Lipscitz number stops decreasing.
% 
% Input:
% * u       ... the system input(column vector)
% * y       ... the system output(column vector)
% * maxlag  ... the max lag space to investigate
% * model   ... optional - 'arx' if we only want to investigate m = l case 
%
% Output:
% * MOIndex ... the maxlag by maxlag matrix containing calculated Lipscitz
%               numbers (Model order index) for each combination of m and l


%% Signature
% Written by Tomaz Sustar
% Based on the algorithm by Xiangdong He and Haruhiko Asada


if(nargin<4), model='unknown'; end;

NN = length(y);                                         % number of samples 
MOIndex = zeros(maxlag);                    % matrix of lipschitz's indexes

for m=1:maxlag,                                 % number of delayed outputs 
  % m,
  for l=0:maxlag,                                 % number of delayed inputs
    
    % if we are investigating arx model srtucture only, we calculate MOIndex only when m = l
    if(strcmp(model, 'arx') && l ~= m), continue; end; 
     
    lag = max(l,m);       % the greater from m, l
    % Because of the lag we can construct only NN - lag input output pairs
    N    = NN-lag;    % number of input - output pairs. 
    p    = floor(0.02*N); % number of Lipschitz quotients used to determine model order index

    [input target] = construct([m l], u, y); % construct regressors and target


    % calculation of Lipschitz quotients
    
     
    Q = zeros(N);      % initialize Q matrix for storing Lipschitz qotients  
    
    for i=1:N-1, 
      % for each input/output pair calculate the their Lipschitz quotients all
      % further inputs/outputs pairs. In this way all possible Lipschitz
      % quotients q(i,j) are calculated.
      
      Q(i,i+1:N)=(target(i)-target(i+1:N)).^2 ./ ...
       sum((repmat(input(i,:), N-i, 1)-input(i+1:N,:)).^2, 2);  
      
    end

    Q_max = Q(Q~=0);                                         % remove zeros
    Q_max = (-sort(-Q_max(:)));        % sort qoutients in descending order
    Q_max = sqrt(Q_max(1:p));                  % take p - largest quotients

    n = m+l;
    MOIndex(m,l+1)=prod(sqrt(n)*Q_max)^(1/p); % calculates order index and stores it to the matrix
    
  end % end for l
  
end % end for m

% draw some figures

if(~strcmp(model, 'arx'))
    
  figure('Name', 'Model order index vs. lag space')
  surf(1:maxlag, 0:maxlag,MOIndex');
  view([-600 40]);
  set(gca, 'Zscale','log');
  set(gca, 'XTick', 1:maxlag)
  set(gca, 'XTick', 1:maxlag)
  xlabel('l - number of past outputs')
  ylabel('m - number of past inputs')
  zlabel('Model Order Index')
end

if(nargin > 4)
  figure(fig);
else
  figure('Name', 'Model order index vs. lag space - arx case')
end
semilogy(diag(MOIndex));
xlabel('m = l - number of past inputs and outputs');
ylabel('Model order index');
set(gca, 'XTick', 1:maxlag);
grid on;




