function [i] = plotgpy(i, t, sys, y, std)
% Function plots the results of the simulation or one-step-ahead
% prediction: target and output together with 95% confidence band.
%
%% Syntax
%  function [i] = plotgpy(i, t, sys, y, std)
%
% Input: 
% * i   ... figure handle, if i==0 function plots in current figure 
% * t   ... time vector (x-axis) 
% * sys ... target output 
% * y   ... predicted output 
% * std ... predicted standard deviation 
% 
% Output: 
% i ... figure handle 
%
% See Also:
% plotg, plotgpe
%
% Examples:
% demo_example_gp_simulation.m
%
%%
% * Written by Dejan Petelin




% check sizes of vectors 
sz = size(t);
if((size(t,1)~=sz(1) | size(sys,1)~=sz(1) |size(y,1)~=sz(1) | size(std,1)~=sz(1))...
        | (size(t,2)~=sz(2) | size(sys,2)~=sz(2) |size(y,2)~=sz(2) | size(std,2)~=sz(2)))
    warning(['figure ', num2str(i), ': vectors: t, tt, y, std must be same size']);
    disp(strcat(['t: ', num2str(size(t))])); 
    disp(strcat(['sys: ', num2str(size(sys))])); 
    disp(strcat(['y: ', num2str(size(y))])); 
    disp(strcat(['std: ', num2str(size(std))])); 
    out = -1; 
    return; 
end

% ix_plot = 1:length(t);   
% % reduce vector if borders are wanted in figure, 
% % e.g. ixplot = 2:length(t)-1;
% 
% xfill = [t(ix_plot); flipdim(t(ix_plot),1)]; 
% yfill = [y(ix_plot)+2*std(ix_plot);flipdim(y(ix_plot)-2*std(ix_plot),1)]; 


% if i==0 use current axis (figure) 
if (i~=0)
figure(i);
end 

%-- dg
iz_cas=t(~isnan(y));
izseki_start=[iz_cas(1);iz_cas(find(diff(find(~isnan(y)))>1)+1)];
izseki_konec=[iz_cas(find(diff(find(~isnan(y)))>1));iz_cas(end)];
% izseki=sum(izseki_start);
% izseki_start=[1;find(diff(find(~isnan(y)))>1)+1];
% izseki_konec=[find(diff(find(~isnan(y)))>1);length(t)];
izseki=length(izseki_start);

plot(t,y, 'k-', 'LineWidth',1); 
hold on 
plot(t,sys, 'k--', 'LineWidth',2); 

for iz=1:izseki
    ix_plot=find(t==izseki_start(iz)) : find(t==izseki_konec(iz));
    %ix_plot=izseki_start(iz):izseki_konec(iz);
    
    xfill = [t(ix_plot); flipdim(t(ix_plot),1)]; 
    yfill = [y(ix_plot)+2*std(ix_plot);flipdim(y(ix_plot)-2*std(ix_plot),1)]; 
    
    fill(xfill, yfill, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
    
end



%-- dg

% fill(xfill, yfill, [7 7 7]/8, 'EdgeColor', [7 7 7]/8);
% hold on 

plot(t,y, 'b-', 'LineWidth',1); 
plot(t,sys, 'g--', 'LineWidth',2); 
%plot(t,y-2*std, 'LineWidth',2, 'Color',[0.7 0.7 0.7]); 
%plot(t,y+2*std, 'LineWidth',2, 'Color',[0.7 0.7 0.7]); 
hold off 
grid on 
xlabel('t'); 
ylabel('Ozone [\mug/m^{3}]'); 
title('GP model simulation')
legend('prediction (\mu)', 'measurement','Location','NorthEast', '\mu \pm 2\sigma'); 
AX=axis; 
AX(1:2)=[t(1) t(end)];  
axis(AX); 





return 






