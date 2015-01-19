%% example_LM_ident
function [lm, Yeq] = example_LM_ident(Ueq, dU, noiseStd, test_flag)

%% Syntax
%  function [lm, Yeq] = example_LM_ident(Ueq, dU, noiseStd, test_flag)

%% Description
% Function which identifies local models for example system (see 
% example.m) in given equilibrium points using Instrumented Variables
% (IV4) algortihm. 
% 
% Inputs:
% Ueq .. equilibrium points where we would like to identify local models 
% dU .. amplitude of the PRBS signal uesd for the perturbation of the
%   system in working points 
% noiseStd .. noise standard deviation of the system's modelled ouput 
% test_flag .. flag for ploting the identification results, 
%   flag=i -> plot figure(i), flag=0 -> no figure 
% Outputs: 
% lm .. structure of local models, each model has the Matlab identification
%   toolbox structure 
% Yeq .. values of the output in corresponding equilibrium points 

%% Examples
% demo_example_gp_data.m

%% See Also
% EXAMPLE, EXAMPLE_DERIVATIVE 


% test_flag tudi figure number 

% constants
Ts = 1;
order = 1; 



    % signals
    startsig = 300;
    multi = 1;
    u = [sig_prbs(7,2)]';
    u = repmat(u,multi,1);
    u = reshape(u,prod(size(u)),1);
    u = dU*u;
    u = Ueq + [zeros(startsig,1); u];
    y = example(u);
    % chopping off one sample
    u = [u(2:end)];
    y = [y(1:end-1)];
    t = 1:length(y);

    % test signal
    uv = [sig_prbs(6,2)]';
    uv = repmat(uv,multi,1);
    uv = reshape(uv,prod(size(uv)),1);
    uv = dU*uv;
    uv = Ueq + [zeros(startsig,1); uv];
    yv = example(uv);
    % chopping off one sample
    %uv = [uv(2:end)];
    uv = [uv(2:end)];
    yv = [yv(1:end-1)];
    tv = 1:length(yv);

    
    % signal processing
    u0 = u(199);
    y0 = y(199);
    
    % noise yes/no 
    if(noiseStd>0)
        noise = normrnd(0, noiseStd, size(y));
        y = y + noise;
    end
    
    u1 = u(startsig+1:end)-u0;
    y1 = y(startsig+1:end)-y0;
    t1 = 1:length(u1);



    % identification
    % lmOrder = [2 2 1];
    lmOrder = [order order 1];
    lm = iv4([y1,u1], lmOrder);
    lm = sett(lm,Ts);
    get(lm) 

    
 Yeq = y0; 




% test
if (test_flag>0)

    % numerical printout - order=1
    dfdy = -lm.a(2);
    dfdu = lm.b(2);

    [dFdYanal dFdUanal Yanal Uanal] = example_derivative(Ueq);
    [dfdy dfdu y0 u0; dFdYanal dFdUanal Yanal Uanal]
    
    %lm
    
    % graphical comparison - identification signal
    
    y2 = sim(lm,u-u0); 
    u2 = u-u0; 
    t2 = Ts:Ts:(length(y2))*Ts; 
    
    y2 = y2 + y0;
    iii = 280:length(y2);
    figure(test_flag*10+1);
    plot(t2(iii),y2(iii),t(iii),y(iii));
    legend('ymdl', 'y', 'Location', 'Best');
    AX = axis(); 
    AX(1:2) = [280 max(t2)]; 
    axis(AX); 
    title(strcat(['Simulation on identification data, U_{eq}=', num2str(u0)])); 
    grid; 
    
    % graphical comparison - validation signal
    yv2 = sim(lm,uv-u0); 
    uv2 = uv-u0; 
    tv2 = Ts:Ts:(length(yv2))*Ts; 
  
    yv2 = yv2 + y0;
    iii = 280:length(yv2);
    figure(test_flag*10+2);
    plot(tv2(iii),yv2(iii),tv(iii),yv(iii));
    legend('ymdl', 'y', 'Location', 'Best');
    AX = axis(); 
    AX(1:2) = [280 max(tv2)]; 
    axis(AX); 
    title(strcat(['Simulation on validation data, U_{eq}=', num2str(u0)])); 
    grid 
    

end


return;



