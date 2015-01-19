function out = sig_prbs(n,m)
% Pseudo-random binary signal (PRBS)
%
%% Syntax
%  out = sig_prbs(n,m)
%
%% Description
% Function generates pseudo-random binary signal (PRBS) with uniform
% amplitude , see R. Iserman: Prozessidentifikation, 1974, Springer. 
%
% Input: 
% * n ... legth of the internal register 
% * m ... number of samples per the width of the narrowest
%         impulse of PRBS  
%
% Output: 
% * out ... PRBS signal values 
% 
% See also:
% sig_prs_minmax 
%
%%

registry1=0:11;
registry2=[0 2 3 4 5 6 7 8 9 10 11 15 20 40 100];
store=ones(n,1);

out(1:m,1)=store(n);

for i=2:2^n-1
    storex=store(registry1(n))+store(registry2(n));
    for j=0:n-2
        store(n-j)=store(n-j-1);
    end
    store(1)=storex;
    if store(1)==2  
        store(1)=0;
    end
    for j=(i-1)*m+1:m*i
        out(j,1)=2*store(n)-1;
    end
end

