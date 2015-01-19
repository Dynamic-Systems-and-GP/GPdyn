%% construct matrix input for training 
function [target,tinput] = construct_train_input(lag,u,y)

%% Syntax
% function [target,input] = construct_train_input(lag,u,y)

%% Description
% Function contructs matrix of regressors for GP model training from input and
% output signals of the system. Shape of the matrix: 
% tinput = [y(1)      ... y(lag)   u(1)       ... u(lag);
%          y(2)       ... y(lag+1) u(2)       ... u(lag+1);
%
%          y(end-lag) ... y(end-1) u(end-lag) ... u(end-1)];
%
% target = [y(lag+1) y(lag+2) ... y(end)]';
%
% Inputs: 
% lag .. order of the system 
% u .. system input signal u, u = [u(1) ... u(end)]' 
% y .. system output signal y, y = [y(1) ... y(end)]' 
% Outputs: 
% target.. target vector for the training routines
% tinput .. matrix of regressors for the training rutines

%% Examples
% demo_example_gp_simulation.m

%% See Also
%
% SIMUL02NAIVE, SIMUL00EXACT, etc. 


if(size(u,1)<size(u,2))
    u = u';
end
if(size(y,1)<size(y,2))
    y = y';
end

len = length(u) - lag + 1;
% input matrix

% y part
for i=1:lag
    tinput(:,i) = y(i:end-lag+i-1);
end
% u part
for i = 1:lag
    tinput(:,lag+i) = u(i:end-lag+i-1);
end
% target

target=y(lag+1:end);

return 
