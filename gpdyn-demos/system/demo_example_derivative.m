function [dFdY, dFdU, Y, U] = demo_example_derivative(U)
% demo_example_derivative - Derives analitically partial derivatives of the
% demo nonlinear system
%% Syntax
%  function [dFdY, dFdU, Y, U] = demo_example_derivative(U)
%
%% Description
% Function which returns partial derivatives df/du and df/dy of the
% nonlinear dynamic system: 
% 
%                y(k)
%   y(k+1) = ---------------  + [u(k)]^3
%            1 + y(k)*y(k)
%
% in equilibrium points. 
% 
% Input:
% U    ... vector, which defines the equilibrium points with the input u.  
%
% Output: 
% dFdY ... partial derivatives df(u,y)/dy
% dFdU ... partial derivatives df(u,y)/du
% Y    ... corresponding f(y) in derivative points
% U    ... corresponding u in derivative points 
%
% See Also
% example.m, demo_example_LM_ident.m
%
% Examples
% demo_example_gp_data.m


if(size(U,1)<size(U,2))
    U = U'; 
end 

for ii=1:length(U)
        u = U(ii); 

        % get y_equilibrium 
        R = roots([1 -u^3 0 -u^3]); 
        for jj=1:3
            if(isreal(R(jj)))
                y = R(jj);       
            end
            break
        end

        % coefficient -a1
        dfdy = (1-y.^2)/(1+y.^2).^2; 
        % coefficient b1 
        dfdu = 3*u.^2; 
        
        Y(ii,1) = y; 
        dFdY(ii,1) = dfdy; 
        dFdU(ii,1) = dfdu; 
        
end 


