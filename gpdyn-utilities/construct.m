function [in, t] = construct(lag, x, y, xind, yind)
% Constructs input matrix and output vector for training or simulation.
%
%% Syntax
% [in,t] = construct(lag, x, y, xind);
% [in,t] = construct([lagy lagx], x, y, xind);
%% Description
% Function contructs matrix of regressors for GP model training from input 
% and output signals of the system. Shape of the matrix:
% 
% in = [y(1)       ... y(lag)   x(1)       ... x(lag);
% 
%       y(2)       ... y(lag+1) x(2)       ... x(lag+1);
%
%       y(end-lag) ... y(end-1) x(end-lag) ... x(end-1)];
%
% If the passed y is empty and yind is not given, the system will be
% considered as FIR. 
% 
% If the passed y is of the size of lag its values will be considered as
% the initial values for simulation (ARX model):
% 
% in = [y(1)       ... y(lag-1)   y(lag)   x(1)       ... x(lag);
% 
%       y(2)       ... y(lag)     0        x(2)       ... x(lag+1);
% 
%       y(lag)     ... 0          0        x(lag)     ... x(lag+lag);
%
%       0          ... 0          0        x(end-lag+1) ... x(end)];
%
% t  = [y(lag+1) y(lag+2) ... y(end)]';
%
% If lag is of a row vector of two integers, first element defines the
% number of past outputs and the secod element defines the nuber of past
% inputs.
%
% Input:
% * lag  ... order of the system
% * x    ... system input signals, x = [x1(1)   ... xN(1);
%                                       x1(2)   ... xN(2);
%                                       x1(end) ... xN(end)],
%            if x is empty, the system will be considered as FIR
% * y    ... system output signal, y = [y(1)    ... y(end)]' or
%            index of system output signal contained in x
% * xind ... indices of signals to be considered (optional)
% * yind ... indices of signals from x to be considered as y
%            output(optional) 
%
% Output:
% * in   ... input matrix for the training rutines
% * t    ... target vector for the training routines
%
% See also:
% simulGPnaive, simulGPmcmc, simulGPexactLIN, simulGPexactSE
%
% Examples:
% demo_example_gp_simulation
%
%%
% Written by D. Petelin, March 2012
% Based on the work of K. Azman

if (nargin < 2)
  error('Error: lag and x are mandatory!');
end
if (nargin < 3)
  y=[];
end

if(isequal (size(lag), [1 2]))
  lagy=lag(1);
  lagx=lag(2);
else
  lagy=lag;
  lagx=lag;
end
maxlag = max(lag);

if (size(x,1) < size(x,2))         % no. of columns bigger than no. of rows
  % no. of regressors > no. of data samples
  r = input('Do your data contain more regressors than data samples? Y/N [N]: ', 's');
  if isempty(r)
    r = 'N';
  end
  if (r == 'N')
    x = x';                                                  % transpose it
  end
end

if (nargin < 5 )
  % y - target
  if (size(y,2) > 1)                                  
    if (size(y,1) > 1)                            % matrix, but expected is
      error('Error: y - number or vector expected!');    % number or vector
    else                                                   % row, therefore
      y = y';                                                % transpose it
    end
  end
                                  % column
  if (length(y) == size(x,1))
    x = [x y];                               % append it to the data matrix
    yind = size(x,2);          % index of target = index of the last column
  elseif (length(y) == lagy)        % just initial values => simulation data
    y = [y; zeros(size(x,1)-lagy+1, 1)];
    x = [x; zeros(1, size(x,2))];
    x = [x y];        
    yind = size(x,2);          % index of target = index of the last column
  elseif (isempty(y))
    yind = -1;                                         % y not given => FIR
  else
    error('Error: Unless y contains inital values, must the lengths of x and y be the same!');   % len. not equal
  end
else   
    %elseif (size(y,1) == 1)
    if (yind >= 1 && yind <= size(x,2))          % index of target is given
      y = x(:,yind);                                   % save target column
    else
      error('Error: yind index out of bounds!');
    end
end

% x - input (indices)
if (nargin > 3 && ~isempty(xind));
  xind = intersect(xind,[1:size(x,2)]);  % consider only reasonable indices
else
  xind = [1:size(x,2)];                             % all available indices
end
if (yind > -1)                                        % if y index is given
  xind = setdiff(xind,yind);                               % remove y index
end

% construct input
in = [];

% y part
if (yind > -1)
  for i = 1 : lagy
    in(:,i) = y(i+maxlag-lagy:end-lagy+i-1);
  end
end
% x part
for i = 1 : lagx
  in = [in x(i+maxlag-lagx:end-lagx+i-1,xind)];
end
% target
if (yind > -1)
  % construct target arx
  t = y(maxlag+1:end);
% else
%   % construct target for time series
%   t = x(lag+1:end);
end


