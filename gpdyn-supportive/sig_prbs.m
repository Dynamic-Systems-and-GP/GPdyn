%% sig_prbs 
function y = sig_prbs(n,m)

%% Syntax
%  function y = sig_prbs(n,m)

%% Description
% Function generates pseudo-random binary signal (PRBS) with uniform
% amplitude, see D. Matko: Identifikacije, 1998, FE Ljubljana, p. 43-45 
% Inputs: 
% n .. legth of the register 
% m .. number of samples per tact of PRBS (the width of the narrowest
% impulse at sampling steps) 
% Output: 
% y .. signal values 


ind1=[0 1 2 3 3 5 4 5 5 7 10];
ind2=[0 2 3 4 5 6 7 8 9 10 11 15 20 40 100];
z=ones(n,1);
for j=1:m
  y(j,1)=z(n);
end
for i=2:2^n-1
  ztemp=z(ind1(n))+z(ind2(n));
  for j=0:n-2
    z(n-j)=z(n-j-1);
  end
  z(1)=ztemp;
  if z(1)==2  
     z(1)=0;
  end
  for j=(i-1)*m+1:m*i
    y(j,1)=2*z(n)-1;
  end
end

