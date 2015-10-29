% demo_example_lmgp_data - demonstration script file for LMGP model identification.  
%
%% Description
% Demo to present how to obtain and compose data for the identification
% of GP model with incorporated local models (LMGP model). 
% 
% There are two options how to obtain local model data: 
% a. analytically (set fIdentifyLm to 0) 
% b. with identification (set fIdentifyLm to 1) 
% 
% It can be used only with the Gaussian covariance function and
% with the white noise model (sum of covSEard and covNoise). 
%
% See Also
% example.m, demo_example_LMGP_training.m, demo_example_LMGP_simulation.m,
% example_LM_ident.m, example_derivative.m
%
%% 
% Written by K. Azman, 2007
%
% Changelog:
%
% 16.2.2015, Martin Stepancic:
%		 	-included the gathering of both equilibrium and non-equilibrium
%			 points in the same matlab script;
%			-changed .mat filenames and their structure contents to: 
%			 {train|valid}_data.{x|u|y}
%
clear;
global flag_LM_data_ident
% LM generation data 
flag_LM_data_ident = 0; % 0=calculated, 1=identified using IV

me = mfilename;                                            % what is my filename
mydir = which(me); mydir = mydir(1:end-2-numel(me));        % where am I located
addpath([mydir 'system']);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 Points at Equilibrium                              %%%%

% local model ident 
ord = 1;    % order of identified LMs 
eq.U = [-1.5 -1.1 -0.7 -0.3 0 0.3 0.7 1.1 1.5]'; 
dU = 0.1;  % depending on diference between the elements of eq.U, valid only if LMs identified!  

noiseStd = 0.1; 
% lm calculated/identified 
if (flag_LM_data_ident)
    % *** local model's parameters identified and not calculated ***    
    
    for ii = 1:length(eq.U)
        [lm{ii},y0] = demo_example_LM_ident(eq.U(ii), dU, noiseStd, 2000+ii);
        
        eq.Y(ii,1) = y0; 
        eq.dfdy(ii,1) = -lm{ii}.a(2); 
        eq.dfdu(ii,1) = -lm{ii}.b(2); 
    end
    eq.lm = lm; 
    
    filename = 'example_data_lmgp_eq_ident';
    
else
    % *** local model's parameters calculated ***
    [eq.dfdy, eq.dfdu, eq.Y] = demo_example_derivative(eq.U); 
    
    filename = 'example_data_lmgp_eq_anal';
   
end

figure('Name', 'EQ Curve');
plot(eq.U,eq.Y);

save(filename,'eq');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                 Non-eqilibrium points                              %%%%
umax = 1.5; 
umin = -umax; 
N_training_points = 200; 
N_validation_points = 200; 
noise_std = 0.05; 

utrain = sig_prs_minmax(N_training_points,2,umin,umax); % create input singal

% simulate example system, to obtain delayed outputs and delayed outputs 
[ytrain, xtrain, utrain] = demo_example(utrain); 
% add noise to training data 
[ytrain, ytrain_no_noise] = add_noise_to_vector(ytrain, noise_std); 
figure('Name', 'Training Input Signal');
t=1:length(utrain);
plot(t, utrain, t, xtrain);
legend('u', 'y');
grid on;

% validation data --> less excited signal 
uvalid = sig_prs_minmax(N_validation_points,5,umin,umax); 
[yvalid, xvalid, uvalid] = demo_example(uvalid); 

figure('Name', 'Validation Input Signal');
plot(uvalid);

% training data 
train_data.u = utrain; 
train_data.x = xtrain; 
train_data.y = ytrain; 

% validation data 
valid_data.u = uvalid; 
valid_data.x = xvalid; 
valid_data.y = yvalid; 

save example_data_lmgp_oeq valid_data train_data umax umin noise_std
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



save(filename,'eq')
rmpath([mydir 'system'])
return
