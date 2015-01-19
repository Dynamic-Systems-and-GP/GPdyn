function [i] = plotgp(i, t, sys, y, std)
% Function plots the results of the simulation or one-step-ahead
% prediction. 
%
%% Syntax
% [i] = plotgp(i, t, sys, y, std);
%
%% Description
% In the upper 2/3 of the figure the output together with 95%
% confidence band and target is plotted, in the lower third of the figure
% the error with corresponding 95% confidence band is plotted. Can be used
% to plot in new or currently active figure. 
%
% Input: 
% * i   ... figure handle, if i==0 function plots in current figure 
% * t   ... time vector (x-axis) 
% * sys ... target output 
% * y   ... predicted output 
% * std ... predicted standard deviation 
% 
% Output: 
% * i ... figure handle 
%
% See Also:
% plotgpy, plotgpe
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

ix_plot = 1:length(t);   
% reduce vector if borders are wanted in figure, 
% e.g. ixplot = 2:length(t)-1;

xfill = [t(ix_plot); flipdim(t(ix_plot),1)]; 
yfill = [y(ix_plot)+2*std(ix_plot);flipdim(y(ix_plot)-2*std(ix_plot),1)]; 


% if i==0 use current axis (figure) 
if (i~=0)
figure(i);
end 

% upper part of figure: y 
subplot(3,1,1:2)
fill(xfill, yfill, [7 7 7]/8, 'EdgeColor', [0 0 0]/8);
hold on 
plot(t,y, 'k-', 'LineWidth',1); 
plot(t,sys, 'k--', 'LineWidth',2); 
%plot(t,y-2*std, 'LineWidth',2, 'Color',[0.7 0.7 0.7]); 
%plot(t,y+2*std, 'LineWidth',2, 'Color',[0.7 0.7 0.7]); 
hold off 
grid on 
xlabel('t'); 
ylabel('y'); 
title('GP model simulation')
legend('\mu \pm 2\sigma', '\mu', 'system','Location','NorthEast'); 
AX=axis; 
AX(1:2)=[t(1) t(end)];  
axis(AX); 

% lower part of figure: e  
subplot(3,1,3)
xfill = [t(ix_plot); flipdim(t(ix_plot),1)]; 
yfill = [2*std(ix_plot);zeros(size(std(ix_plot)))]; 
fill(xfill, yfill, [7 7 7]/8, 'EdgeColor', [0 0 0]/8);
% plot(t,2*std, 'LineWidth',2, 'Color',[0.7 0.7 0.7]); 
hold on 
plot(t,abs(y-sys), 'k-', 'LineWidth',1); 
hold off 
legend('2\sigma', '|e|', 'Location','NorthEast'); 
grid 
xlabel('t'); 
ylabel('e'); 
AX=axis; 
AX(1:2)=[t(1) t(end)]; 
axis(AX); 



return 






