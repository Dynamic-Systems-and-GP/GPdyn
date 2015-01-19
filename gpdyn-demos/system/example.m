%% example 
function [Y, X, U] = example(U, varargin)

%% Syntax
%  function [Y, X, U] = example(U, varargin)

%% Description
% Function simulating the nonlinear dynamic system: 
% 
%                y(k)
%   y(k+1) = ---------------  + [u(k)]^3
%            1 + y(k)*y(k)
%
% which is used to demonstrate the use of this toolbox. 
% [1] K.S. Narendra and K. Parthasarathy. Identification
% and Control of Dynamical Systems Using Neural Networks, 
% IEEE Transactions on NN, Vol.1 No. 1, 4-27, 1990.
%
% Usage: 
% (1) [Y, X, U] = example(u, y0), u is vector, y0 is optional 
% (2) [Y, X, U] = example(U, y0), U is Nx2 matrix 
% Inputs:
% u .. vector of inputs 
% U .. matrix U=[x u], where x is vector of past states nad u is vector of
%   inputs
% y0 .. intial state of the system, optional 
% Outputs: 
% Y .. vector of outputs 
% X .. vector of past outputs: X(k) = Y(k-1)
% U .. vector of inputs 
% 
% 

%% Examples
% demo_example_gp_data.m

%% See Also
% EXAMPLE_DERIVATIVE, EXAMPLE_LM_IDENT


if(nargin==1)
    y0 = 0; 
else 
    y0 = varargin{1}; 
end 

if (size(U,2)>2 | size(U,2)<1)
    error('input U must have size 1 or 2'); 
elseif (size(U,2)==2)   
    % first coloumn = y
    % second coloumn = u 
    X = U(:,1); 
    U = U(:,2); 
    Y = X./(1+X.^2) + U.^3; 
else % size(U,2)=1 
    % step1 
    X = y0; 
    % steps 
    for ii=1:length(U)
        Y(ii) = X(ii)/(1+X(ii).^2) + U(ii).^3; 
        X(ii+1) = Y(ii); 
    end 
    
    X = X(1:end-1); 
    X = X'; Y = Y'; 
end 
    
    

