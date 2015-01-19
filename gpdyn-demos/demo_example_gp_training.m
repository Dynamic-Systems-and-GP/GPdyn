%% demo_example_gp_training

%% Description
% Demo to present the how to train (=identify) the GP model, which
% describes the nonlinear dynamic system. 

%% See Also
% EXAMPLE, DEMO_EXAMPLE_GP_DATA, DEMO_EXAMPLE_GP_SIMULATION, TRAINGP, GPR 
 
clear;
close all;

load example_data 

yeq = example_train_data.yeq; 
ueq = example_train_data.ueq; 
utrain = example_train_data.utrain; 
xtrain = example_train_data.xtrain; 
ytrain = example_train_data.ytrain; 

% build training data (delayed outputs y first) 
input = [yeq ueq; xtrain utrain]; 
target = [yeq; ytrain]; 

% covariance function: SE + white noise 
covfunc = {'covSum',{'covSEard','covNoise'}}; 

% to get the number of hyperparameters you could use: feval(covfunc{:})
D = size(input,2); 
logtheta0 = -ones(D+2,1); 


% training 
[logtheta, flogtheta, i] = traingp(covfunc,input, target, logtheta0);

% validation on identification data 
[ytest S2test] = gpr(logtheta, covfunc, input, target, input);

figure 
plot(input(:,1),target,'b*'); 
hold on 
plot(input(:,1),ytest,'ko'); 
plot(input(:,1),ytest-2*sqrt(S2test),'k.'); 
plot(input(:,1),ytest+2*sqrt(S2test),'k.'); 
hold off 
title('Validation on training data')
xlabel('y(k-1)')
ylabel('y(k)')
legend('target','ygp','ygp\pm2\sigma') 

save example_trained 



return; 



