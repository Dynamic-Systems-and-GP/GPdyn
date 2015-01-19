%% detect_fault 
function [index_bool, index_val, areas] = detect_fault(residuals, var_residuals, NW, c_max)

%% Syntax
%  function [index_bool, index_val] = detect_fault(residuals, var_residuals, NW, c_max)

%% Description
% Function for detecting the fault, where the system is modelled with the
% GP model. It takes into account predicted mean as well as the predicted
% variance. Greater window size can be selected to increase robustness on
% the account of lower sensibility. 
% Note: Function is not compatible with the current version of the GPML 
% toolbox. 
% [1] Dj.Juricic, J.Kocijan,Fault detection based on Gaussian process model,
% I.Troch, F.Breitenecker (eds.), Proceedings of the 5th Vienna
% Symposium on Mathematical Modelling – MathMod, Wien, 2006
% 
% Inputs: 
% residuals .. [eps(1) .. eps(n)], where eps(k) = y(k)-yhat(k) = y(k)-m(k)
% var_residuals .. [v(1) .. v(n)], where w(k) = var(eps(k))  predicted
%       variance associated with yhat(k), i.e. GP prediction at time step k
% NW .. size of window, used for calculating the index and deciding about fault
% c_max .. maximal significance level still reflecting in deciding for H0
%   (no fault)
% Outputs: 
% index_bool: 1 for H1 (fault) and 0 for H0 (no fault)
% index_val: significance levels, corresponding to index_bool  


%% See Also
% DETECT_FAULT_VALIDITY, SIMUL00MCMC



% flag=1: while debugging and testing
flag_plot = 0;

if (size(residuals)~=size(var_residuals))
    error('vectors residuals and var_residuals must be of equal length');
elseif (NW>length(residuals))
    error('window size NW must be smaller or equal then the length of residuals');
end


% classificator 
for ii=NW:length(residuals)

    SIGMA = diag(var_residuals(ii:-1:ii-NW+1));
    res = residuals(ii:-1:ii-NW+1);
    % from paper [1]  
    mi(ii,1) = abs(ones(1,NW)*inv(SIGMA)*res)/(ones(1,NW)*inv(SIGMA)*ones(NW,1));
    % mi2 - from paper [1]
    c(ii,1) = abs(ones(1,NW)*inv(SIGMA)*res)/sqrt((ones(1,NW)*inv(SIGMA)*ones(NW,1)));
end


t = [0:length(residuals)-1]';
if (flag_plot==1)
    figure(550)
    subplot(3,1,1)
    plot(t,residuals)
    legend('residuals',2)
    grid 
    axis([t(1) t(end) floor(min(residuals)) ceil(max(residuals))]);
    subplot(3,1,2)
    plot(t,var_residuals)
    legend('var',2) 
    grid 
    subplot(3,1,3)
    plot(t, mi)
    legend('mi',2)
    grid 
    
    
    figure(551) 
    plot(t,c,t,c_max*ones(size(t)))
    title('value of cost c ')
    grid 
end


index_val = c; 
index_bool = ge(c, c_max*ones(size(mi))); 

index_bool_start = []; 
index_bool_stop = []; 
flag = 0; 
jj = 1; 
for ii=1:length(t)
    if((flag==0) & (index_bool(ii)==1))
        index_bool_start(jj) = ii; 
        flag = 1; 
    end 
    if((flag==1) & (index_bool(ii)==0))
        index_bool_stop(jj) = ii-1; 
        flag = 0; 
        jj = jj+1; 
    end 
end 
areas.start = index_bool_start; 
areas.stop = index_bool_stop; 

    




return;

