%% construct matrix input for simulation 
function xt = construct_simul_input(lag,u,x0)

%% Syntax
% function xt = construct_simul_input(lag,u,x0)

%% Description
% Function contructs matrix for GP model simulation from input signal and
% starting values of output of the system. Shape of the matrix: 
% xt = [y(1)... y(lag) u(1)      ... u(lag);
%       y(2)... 0      u(2)      ... u(lag+1);
%
%       0   ... 0      u(end-lag)... u(end-1)];
%
% Inputs: 
% lag .. order of the system 
% u .. system input signal u, u = [u(1) ... u(end)]' 
% x0 .. starting values of system state (=output),  x0 = [y(1) ... y(lag)] 
% Outputs: 
% xt .. matrix used as the input to simulation rutines

%% Examples
% demo_example_gp_simulation.m

%% See Also
%
% SIMUL02NAIVE, SIMUL00EXACT, etc. 


if(size(u,1)<size(u,2))
    u = u';
end

len = length(u) - lag + 1;
% matrix
xt = zeros(len, 2*lag);
% y part
for i=1:lag
    xt(i,1:lag-i+1) = x0(i:end);
end
% u part
for i = 1:lag
    xt(:,lag+i) = u(i:end-lag+i);
end

return 
