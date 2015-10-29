function [in, t] = construct(lag, u, y, uind, yind)
% Constructs the input matrix and the output vector for training or simulation.
%
%% Syntax
% [in,t] = construct(lag, u, y, uind);
% [in,t] = construct([lagy lagu], u, y, uind);
%% Description
% Function contructs matrix of regressors for GP model training from input 
% and output signals of the system. Shape of the matrix:
% 
% in = [y(1)       ... y(lag)   u(1)       ... u(lag);
% 
%       y(2)       ... y(lag+1) u(2)       ... u(lag+1);
%
%       y(end-lag) ... y(end-1) u(end-lag) ... u(end-1)];
%
% If the passed y is empty and yind is not given, the system will be
% considered as FIR. 
% 
% If the passed y is of the size of lag its values will be considered as
% the initial values for simulation (ARX model):
% 
% in = [y(1)       ... y(lag-1)   y(lag)   u(1)       ... u(lag);
% 
%       y(2)       ... y(lag)     0        u(2)       ... u(lag+1);
% 
%       y(lag)     ... 0          0        u(lag)     ... u(lag+lag);
%
%       0          ... 0          0        u(end-lag+1) ... u(end)];
%
% t  = [y(lag+1) y(lag+2) ... y(end)]';
%
% If lag is of a row vector of two integers, first element defines the
% number of past outputs and the secod element defines the nuber of past
% inputs.
%
% Input:
% * lag  ... the order of the system
% * u    ... system input signals, u = [u1(1)   ... uN(1);
%                                       u1(2)   ... uN(2);
%                                       u1(end) ... uN(end)],
% * y    ... the system output signal, y = [y(1)    ... y(end)]' or
%            index of system output signal contained in u
% * uind ... indices of signals to be considered (optional)
% * yind ... indices of signals from u to be considered as y
%            output(optional) 
%
% Output:
% * in   ... the input matrix for the training rutines
% * t    ... the target vector for the training routines
%
% See also:
% simulGPnaive, simulGPmc, simulGPexactLIN, simulGPexactSE
%
% Examples:
% demo_example_gp_simulation
%
%%
% Written by D. Petelin, March 2012
% Based on the work of K. Azman

if (nargin < 2)
  error('Error: lag and u are mandatory!');
end
if (nargin < 3)
  y=[];
end

if(isequal (size(lag), [1 2]))
  lagy=lag(1);
  lagu=lag(2);
else
  lagy=lag;
  lagu=lag;
end
maxlag = max(lag);

if (size(u,1) < size(u,2))         % no. of columns bigger than no. of rows
  % no. of regressors > no. of data samples
  r = input('Do your data contain more regressors than data samples? Y/N [N]: ', 's');
  if isempty(r)
    r = 'N';
  end
  if (r == 'N')
    u = u';                                                  % transpose it
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
  if (length(y) == size(u,1))
    u = [u y];                               % append it to the data matrix
    yind = size(u,2);          % index of target = index of the last column
  elseif (length(y) == lagy)        % just initial values => simulation data
    y = [y; zeros(size(u,1)-lagy+1, 1)];
    u = [u; zeros(1, size(u,2))];
    u = [u y];        
    yind = size(u,2);          % index of target = index of the last column
  elseif (isempty(y))
    yind = -1;                                         % y not given => FIR
  else
    error('Error: Unless y contains inital values, must the lengths of u and y be the same!');   % len. not equal
  end
else   
    %elseif (size(y,1) == 1)
    if (yind >= 1 && yind <= size(u,2))          % index of target is given
      y = u(:,yind);                                   % save target column
    else
      error('Error: yind index out of bounds!');
    end
end

% u - input (indices)
if (nargin > 3 && ~isempty(uind));
  uind = intersect(uind,[1:size(u,2)]);  % consider only reasonable indices
else
  uind = [1:size(u,2)];                             % all available indices
end
if (yind > -1)                                        % if y index is given
  uind = setdiff(uind,yind);                               % remove y index
end

% construct input
in = [];

% y part
if (yind > -1)
  for i = 1 : lagy
    in(:,i) = y(i+maxlag-lagy:end-lagy+i-1);
  end
end
% u part
for i = 1 : lagu
  in = [in u(i+maxlag-lagu:end-lagu+i-1,uind)];
end
% target
if (yind > -1)
  % construct target arx
  t = y(maxlag+1:end);
end


