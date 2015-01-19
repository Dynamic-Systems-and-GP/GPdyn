function [X, fX, i] = minimizeDE(X, f, itermax, varargin)
% Minimize a multivariate function using differential evolution.
% 
%% Syntax
% [X, fX, i] = trainDEgp(X, f, itermax, P1, P2, P3, ... )
%
%% Description
% Minimization of a user-supplied function using the differential with
% respect to X(1:dim), evolution (DE) algorithm. DE works best if 
% [minbound, maxbound] covers the region where the global minimum is 
% expected. DE is also somewhat sensitive to the choice of the stepsize F. 
% A good initial guess is to choose F from interval [0, 2], e.g. 1.0. CR,
% the crossover probability constant from interval [0, 1] helps to maintain
% the diversity of the population. Only separable problems do better with 
% CR close to 0. If the parameters are correlated, high values of CR work 
% better. The reverse is true for no correlation.
%
% The number of population members npop is also not very critical. A good 
% initial guess is 10*dim. Depending on the difficulty of the problem npop 
% can be lower than 10*dim or must be higher than 10*dim to achieve 
% convergence.
%
% minimize_de is a vectorized variant of DE which, however, has a property 
% which differs from the original version of DE: the random selection of 
% vectors is performed by shuffling the population array. Hence a certain 
% vector can't be chosen twice in the same term of the perturbation 
% expression. Due to the vectorized expressions minimize_de executes fairly
% fast in MATLAB's interpreter environment.
%
% Input:
% * X       ... initial guess; may be of any type, including struct and cell 
%                array
% * f       ... the name or pointer to the function to be minimized. The 
%               function f must return the value of the function.
% * itermax ... number of iterations.
% * P1, P2, ... parameters are passed to the function f.
%
% Output:
% * X       ... the returned solution
% * fX      ... vector of function values indicating progress made
% * i       ... number of iterations (generations) used at termination.
%
% Examples:
%
% See also:
%
%%
%
% Author: Dejan Petelin, based on an algorithm by Kenneth Price and Rainer
% Storm 
% 
% Note:
% This program is free software; you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation; either version 1, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. A copy of the GNU 
% General Public License can be obtained from the 
% Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.



if (itermax <= 0) % validation
  itermax = 1000;
  fprintf('Number of iterations (length) should be > 0; ');
  fprintf('set to default value 200.\n');
end

% parameters
npop = 50;                                  % Number of population memebers
F = 0.5;                               % DE-stepsize F from interval [0, 2]
CR = 0.5;          % crossover probability constant CR from interval [0, 1]
dim = length(unwrap(X));       % number of params of the objective function
minbound = repmat(-8, 1, dim);     % vector of lower bounds of initial pop.
maxbound = repmat(7, 1, dim);      % vector of upper bounds of initial pop.
bounds = 1;                  % vector of lower bounds of initial population
strategy = 1;                                         % strategy: DE/rand/1
refresh = 1;                               % print results only on progress

% validation
if (npop < 5)
  npop = 5;
  fprintf('Number of population (npop) increased to minimal value 5\n');
end
if ((F < 0) || (F > 1))
  F = 0.5;
  fprintf('F should be from interval [0,2]; set to default value 0.5\n');
end
if ((CR < 0) || (CR > 1))
  CR = 0.5;
  fprintf('CR should be from interval [0,1]; set to default value 0.5\n');
end

% intialization
pop = zeros(npop, dim);                      % initialize pop to gain speed
pop(1, :) = unwrap(X);                              % passed initial values
for k = 2 : npop   % random values between min and max values of the params
  pop(k, :) = minbound + rand(1, dim) .* (maxbound - minbound);
end

popold = zeros(size(pop));                              % toggle population
membest = zeros(1, dim);                      % best population member ever
membestit = zeros(1, dim);            % best population member in iteration
nfeval = 0;                                % number of function evaluations

% evaluation
ibest = 1;                             % start with first population member
val(1) = feval(f, rewrap(X, pop(ibest, :)), varargin{:});
valbest = val(1);                    % best objective function value so far
fX = valbest;
nfeval = nfeval + 1;
for k = 2 : npop                              % check the remaining members
  val(k) = feval(f, rewrap(X, pop(k, :)), varargin{:});
  nfeval = nfeval + 1;
  if (val(k) < valbest)
     ibest = k;                                         % save its location
     valbest = val(k);
  end   
end
membestit = pop(ibest, :);               % best member of current iteration
valbestit = valbest;                      % best value of current iteration
membest = membestit;                                     % best member ever
fX = [fX' valbest]';

pm1 = zeros(npop, dim);                    % initialize population matrix 1
pm2 = zeros(npop, dim);                    % initialize population matrix 2
pm3 = zeros(npop, dim);                    % initialize population matrix 3
pm4 = zeros(npop, dim);                    % initialize population matrix 4
pm5 = zeros(npop, dim);                    % initialize population matrix 5
bm = zeros(npop, dim);                         % initialize membest  matrix
ui = zeros(npop, dim);       % intermediate population of perturbed vectors
mui = zeros(npop, dim);                  % mask for intermediate population
mpo = zeros(npop, dim);                           % mask for old population
rot = 0 : 1 : npop - 1;                  % rotating index array (size npop)
rotd = 0 : 1 : dim - 1;                   % rotating index array (size dim)
rt = zeros(npop);                            % another rotating index array
rtd = zeros(npop);                 % rotating index array for exp crossover
a1 = zeros(npop);                                             % index array
a2 = zeros(npop);                                             % index array
a3 = zeros(npop);                                             % index array
a4 = zeros(npop);                                             % index array
a5 = zeros(npop);                                             % index array
ind = zeros(4);
meanv = ones(npop, dim);

% DE minimization
iter = 1;
while (iter < itermax)
  prog = 0;
  popold = pop;                                   % save the old population
  ind = randperm(4);                                  % index pointer array

  a1 = randperm(npop);                       % shuffle locations of vectors
  rt = rem(rot + ind(1), npop);        % rotate indices by ind(1) positions
  a2 = a1(rt + 1);                                % rotate vector locations
  rt = rem(rot + ind(2), npop);
  a3 = a2(rt + 1);
  rt = rem(rot + ind(3), npop);
  a4 = a3(rt + 1);
  rt = rem(rot + ind(4), npop);
  a5 = a4(rt + 1);

  pm1 = popold(a1, :);                              % shuffled population 1
  pm2 = popold(a2, :);                              % shuffled population 2
  pm3 = popold(a3, :);                              % shuffled population 3
  pm4 = popold(a4, :);                              % shuffled population 4
  pm5 = popold(a5, :);                              % shuffled population 5

  mui = rand(npop, dim) < CR;  % all random numbers < CR are 1, 0 otherwise
  
  for k = 1 : npop                 % population filled with the best member
    bm(k, :) = membestit;                           % of the last iteration
  end
  
  % Insert this if you want exponential crossover
  %mui = sort(mui)';	              % transpose, collect 1's in each column
  %for k = 1 : npop
  %  n = floor(rand * dim);
  %  if (n > 0)
  %     rtd = rem(rotd + n, dim);
  %     mui(:,k) = mui(rtd + 1, k);                   %rotate column k by n
  %  end
  %end
  %mui = mui';			                                       % transpose back
  % End of exponential crossover
  
  mpo = mui < 0.5;                                    % inverse mask to mui
  
  % strategies
  if (strategy == 1)                                            % DE/rand/1
    ui = pm3 + F * (pm1 - pm2);                    % differential variation
    ui = popold .* mpo + ui .* mui;                             % crossover
    origin = pm3;
  elseif (strategy == 2)                               % DE/local-to-best/1
    ui = popold + F * (bm - popold) + F * (pm1 - pm2);
    ui = popold .* mpo + ui .* mui;
    origin = popold;
  elseif (strategy == 3)                            % DE/best/1 with jitter
    ui = bm + (pm1 - pm2) .* ((1-0.9999) * rand(npop, dim) + F);               
    ui = popold .* mpo + ui .* mui;
    origin = bm;
  elseif (strategy == 4)                 % DE/rand/1 with per-vector-dither
    f1 = ((1 - F) * rand(npop, 1) + F);
    for k = 1 : dim
      pm5(:, k) = f1;
    end
    ui = pm3 + (pm1 - pm2) .* pm5;                 % differential variation
    origin = pm3;
    ui = popold .* mpo + ui .* mui;                             % crossover
  elseif (strategy == 5)                 % DE/rand/1 with per-vector-dither
    f1 = ((1 - F) * rand + F);
    ui = pm3 + (pm1 - pm2) * f1;                   % differential variation
    origin = pm3;
    ui = popold .* mpo + ui .* mui;                             % crossover
  else                                                % either-or-algorithm
    if (rand < 0.5);                                            % Pmu = 0.5
      ui = pm3 + F * (pm1 - pm2);                  % differential variation
      origin = pm3;
    else                                       % use F-K-Rule: K = 0.5(F+1)
      ui = pm3 + 0.5 * (F + 1.0) * (pm1 + pm2 - 2 * pm3);
    end
    ui = popold .* mpo + ui .* mui;                             % crossover     
  end
  
  % select which vectors are allowed to enter the new population
  for k = 1 : npop
    if (bounds)                        % if boundary constraints are needed
      for j = 1 : dim                % boundary constraints via bounce back
        if (ui(k, j) > maxbound(j))
          ui(k, j) = maxbound(j) + rand * (origin(k, j) - maxbound(j));
        end
        if (ui(k, j) < minbound(j))
          ui(k, j) = minbound(j) + rand * (origin(k, j) - minbound(j));
        end   
      end
    end
  
    % check cost of competitor
    try
      valtemp = feval(f, rewrap(X, ui(k, :)), varargin{:});
    catch
      valtemp = Inf;
    end
    nfeval = nfeval + 1;
    if (valtemp < val(k))
      pop(k, :) = ui(k, :);               % replace old vector with new one
      val(k) = valtemp;                        % save value in "cost array"
      if (valtemp < valbest)       % update valbest only in case of success
        prog = 1;
        valbest = valtemp;                                 % new best value
        membest = ui(k, :);                % new best parameter vector ever
        X = rewrap(X, membest);
      end
    end
  end

  if (prog)
    fX = [fX' valbest]';
  end
  
  membestit = membest;                      % best member of this iteration

  % output
  if (refresh == 1)
    if (prog || (iter == 1))
      fprintf('Iteration %6i;  Value %4.6e\n', iter, valbest);
    end
  elseif (refresh == 2)
    fprintf('Iteration %6i;  Value %4.6e\n', iter, valbest);
  end

  iter = iter + 1;
end
i = iter;

% unwrap
% Extract the numerical values from "s" into the column vector "v". The
% variable "s" can be of any type, including struct and cell array.
% Non-numerical elements are ignored. See also the reverse rewrap.m. 
function v = unwrap(s)

v = [];   
if isnumeric(s)
  v = s(:);                    % numeric values are recast to column vector
elseif isstruct(s)
  v = unwrap(struct2cell(orderfields(s)));   % alphabetize, convert to cell
elseif iscell(s)
  for i = 1:numel(s)         % cell array elements are handled sequentially
    v = [v; unwrap(s{i})];
  end
end                                               % other types are ignored

% rewrap
% Map the numerical elements in the vector "v" onto the variables "s" which
% can be of any type. The number of numerical elements must match; on exit 
% "v" should be empty. Non-numerical entries are just copied. See also
% unwrap.m.
function [s v] = rewrap(s, v)

if isnumeric(s)
  if numel(v) < numel(s)
    error('The vector for conversion contains too few elements')
  end
  s = reshape(v(1:numel(s)), size(s));        % numeric values are reshaped
  v = v(numel(s)+1:end);                    % remaining arguments passed on
elseif isstruct(s) 
  [s p] = orderfields(s); p(p) = 1:numel(p);  % alphabetize, store ordering  
  [t v] = rewrap(struct2cell(s), v);             % convert to cell, recurse
  s = orderfields(cell2struct(t,fieldnames(s),1),p);    % convert to struct
elseif iscell(s)
  for i = 1:numel(s)         % cell array elements are handled sequentially 
    [s{i} v] = rewrap(s{i}, v);
  end
end                                         % other types are not processed