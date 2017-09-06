%% demo_example_gp_training

%% Description
% This demo presents how to train, i.e., identify, GP model, which
% describes the nonlinear dynamic system.

%% See Also
% demo_example, demo_example_gp_data, demo_example_gp_simulation,
% demo_example_norm, gp, gpx, gp_initial, trainGParx, trainGPoe 
 
clear all;
close all;

% load data from file
load example_data 

% Build training data (delayed outputs y first) 
input = [xtrain utrain]; 
target = [ytrain]; 

% Define covariance function: the squared exponential covariance function with
% ARD
cov = @covSEard; 

% Define covariance function: the Gaussian likelihood function.
lik = @likGauss;

% Define mean function: the zero mean function.
mean = @meanZero;

% Define inference method: the exact Inference
inf= @infExact;

%% Setting initial hyperparameters
% For hyperparameters a structure array is used. The structure array has to be of the
% following shape:
% hyp = 
%      cov: [] - covariance function parameters
%      lik: [] - likelyhood function parameters
%      mean: [] -  mean function parameters
% 
% To get the number of hyperparameters you may
% use: eval_func(cov),  eval_func(cov) or eval_func(mean).
%

D = size(input,2); % Input space dimension
hyp0.cov  = -ones(D+1,1); 

% Define the likelihood hyperparameter. In our case this parameter is the noise
% parameter.
hyp0.lik=log(0.1);
 
hyp0.mean = []; 

%% gp_initial
% We can also use the function gp_initial to find initial hperparameters.
% This function returns the best set of n random sets of hyperparameter values. 
% As score it uses a log marginal likelihood. The number of parameters
% is adjusted to the current covariance, likelihood and mean function.% 

% Set between which bounds the best set of hyperparameters will be estimated.
bounds=[-7,8];

% Find initial hyperparameters:
hyp0_lin = gp_initial(bounds, inf, mean, @covLINard, lik, input, target);

% For further use we will train another GP model with linear covariance
% function(@covLINard).


%% Training 
disp('Identification of GP model');
[hyp, flogtheta, i] = trainGParx(hyp0, inf, mean, cov, lik, input, target);

disp('Training using Differential Evolution minimization algorithm with default value of iterations:');
[hyp_lin, flogtheta_lin, i] = trainGParx(hyp0_lin, inf, mean, @covLINard, lik, input, target, @minimizeDE);

disp('Training using Output Error algorithm');
[hyp_oe, flogtheta, i,simy,simse2] = trainGPoe(hyp0, inf, mean, cov, lik, input, target);

figure('Name','Iterations of simulated response from Output Error algorithm.');
subplot(2,1,1);hold on;
plot(simy(:,1),'r');                
cl=linspace(1,0,size(simy,2));
for i=2:5:size(simy,2)
	plot(simy(:,i),'color',[0,0,cl(i)]);
end
ylabel('response');legend({'measured response'});
subplot(2,1,2);hold on;
for i=2:5:size(simy,2)
	plot(1:size(simse2,1),2*sqrt(simse2(:,i)),'color',[0,cl(i)/2,cl(i)/2]);
end
ylabel('2 std. devs.');


%% Validation (Regression)
disp('validation on identification data');
[ytest S2test] = gp(hyp, inf, mean, cov, lik, input, target, input);

% plot
t = [0:length(input)-1]';
f1=figure('Name', 'Validation on Identification Data');
plotgp(f1,t,target, ytest, sqrt(S2test));

disp('validation on validation dataset (regression)');
[ytest2 S2test2] = gp(hyp, inf, mean, cov, lik, input, target, [xvalid uvalid]);

%plot
t = [0:length(uvalid)-1]';
f2=figure('Name', 'Validation on Validation Data (Regression)');
plotgp(f2,t,yvalid, ytest2, sqrt(S2test2));

% save trained GP model to file 
save example_trained hyp hyp_lin hyp_oe inf mean cov lik input target



return; 



