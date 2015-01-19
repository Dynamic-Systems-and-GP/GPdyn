%% simul00naive_hystOut
function [y, s2, hyst] = simul00naive_hystOut(logtheta, covfunc, input, target, xt,lag, hystProp) 

%% Syntax
% function [y, s2, hyst] = simul00naive_hystOut(logtheta, covfunc, input,
%   target, xt,lag, hystProp) 

%% Description
% Template for function for "naive" (i.e. without propagation of variance)
% simulation of the GP model with hysteresis on the output. Uses routine
% gpr and works similar to SIMUL00NAIVE. 
% See K. Azman, J. Kocijan, Identifikacija dinamiènega sistema s histerezo
% z modelom na osnovi Gaussovih procesov. In: B. Zajc, A. Trost (eds.) 
% Zbornik štirinajste mednarodne Elektrotehniške in raèunalniške konference
% ERK 2005, 26.-28.09., Portorož, Slovenija, 2005, Vol. A, pp. 253--256. 
% (in Slovene)
% 
% Inputs: 
% loghteta .. optimized hyperparameters 
% covfunc .. specified covariance function, see help covFun for more info 
% input .. input part of the training data,  NxD matrix, D-th regressor
%   represents the state of hysteresis 
% target .. output part of the training data (ie. target), Nx1 vector 
% xt .. input matrix for simulation, kxD vector, see
%   construct_input_matrix.m for more info 
% lag .. the order of the model (number of used lagged outputs) 
% hystProp .. hysteresis properties: 
%   .limits .. [xmin xmax] values, where the state of hysteris changes 
%   .DYH .. change of the output due to hysteresis 
%   .state .. initial state of hysteresis 
% Outputs: 
% y .. mean predicted output 
% s2 .. corresponding variances with added white noise v0 
% hyst .. corresponding states of the hysteresis 

%% See Also
% GPR, SIMUL00NAIVE


fun_name = 'simul00naive_hystOut'; 


hyst_low = hystProp.limits(1);
hyst_high = hystProp.limits(2);
hystState = xt(1,end);

[NN,D] = size(xt);

% first step 
test = xt(1,:);
hyst(1) = hystState;
[y(1), s2(1)] = gpr(logtheta, covfunc, input, target, test);


for k=2:NN
    
    if (mod(k,50) == 0)
        disp([fun_name, ', step: ', int2str(k), '/', int2str(NN)]);
    end    
    
    % hysteresis state     
    if(hystState==-1 & mu(k-1)>(hyst_high-DYH))          
        hystState = 1;         
    elseif(hystState==1 & mu(k-1)<(hyst_low+DYH))        
        hystState = -1; 
    end
    hyst(k) = hystState; 
    
    
    if (k>lag)
        test = [y(k-lag:k-1) xt(k, lag+1:end)];
    elseif (k<=lag)
        test = [xt(k, 1:lag-k+1) y(1:k-1) xt(k, lag+1:end)];
    end
    test(D) = [hystState];
    
    
    [y(k), s2(k)] = gpr(logtheta, covfunc, input, target, test);
    
end

% adding noise variance v0
s2 = s2 + exp(logtheta(end))^2;


% transform into coloumn vectors 
y = y'; 
s2 = s2'; 
hyst = hyst'; 

