%% demo_example_gp_data 

%% Description
% Demo to present how to obtain and compose data for the identification
% of the GP model. 

%% See Also
% EXAMPLE, DEMO_EXAMPLE_PRESENT, DEMO_EXAMPLE_GP_TRAINING, 
% DEMO_EXAMPLE_GP_SIMULATION 


clear;
close all;

% data parameters 
umax = 1.5; 
umin = -umax; 
N_training_points = 100; 
N_validation_points = 200; 
noise_std = 0.05; 


% training data --> eq curve AND samples very excited signal 
ueq = [-1.5:0.5:1.5]';
for ii=1:length(ueq)
    ueq0 = ueq(ii);
    ytmp = example(repmat(ueq0,400,1));
    yeq(ii,1) = ytmp(end);      
end 

utrain = sig_prs_minmax(N_training_points,2,umin,umax); 
[ytrain, xtrain, utrain] = example(utrain); 
% add noise to training data 
[ytrain, ytrain_no_noise] = add_noise_to_vector(ytrain, noise_std); 

% validation data --> less excited signal 
uvalid = sig_prs_minmax(N_validation_points,1,umin,umax); 
[yvalid, xvalid, uvalid] = example(uvalid); 

% training data 
example_train_data.ueq = ueq; 
example_train_data.yeq = yeq; 
example_train_data.utrain = utrain; 
example_train_data.xtrain = xtrain; 
example_train_data.ytrain = ytrain; 

% validation data 
example_valid_data.uvalid = uvalid; 
example_valid_data.xvalid = xvalid; 
example_valid_data.yvalid = yvalid; 

save example_data example_valid_data example_train_data umax umin noise_std 




return; 



