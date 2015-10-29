%% demo_example_gp_norm

%% Description
% If different inputs are measured on different scales a normalisation can
% be applied to achieve that all inputs have equal impact on the output.
% There are three functions available in this toolbox to assist with
% normalisation. This demo shows their usage.

%% See Also
% preNorm, postNorm, postNormVar, gp_initial, construct,
% demo_example_gp_training, demo_example_gp_simulation

load example_data                                          % load test data                                         

%% Training

% Build training data (delayed outputs y first) 
input = [xtrain utrain]; 
target = [ytrain]; 

% set input functions
cov = @covSEard; lik = @likGauss; mean = @meanZero; inf = @infExact;

% normalize training data
[nInput,inputMin,inputMax, nTarget, targetMin, targetMax] = preNorm(input, target);                         

% find initial hyperparameters using the normalized data
hyp0n = gp_initial([-1, 1], inf, mean, cov, lik, nInput, nTarget);

% train model
[nHyp, flogtheta, i] = trainGParx(hyp0n, inf, mean, cov, lik, nInput, nTarget);

%% Simulation

% prepare input matrix for simulation
lag = 1;
y0 = 0;
test = construct(lag, uvalid, y0);

% normalize test data using minimums and maximums computed previously 
nTest = preNorm(test,inputMin,inputMax);

% exact simulation using the normalized test data and model
[nY, nS2] = simulGPexactSE(nHyp, inf, mean, cov, lik, nInput, nTarget, nTest, lag);

%% Representing results

% unnormalize simulation mean output
y = postNorm(nY, targetMin, targetMax);

% unnormalize variance
s2 = postNormVar(nS2, targetMin, targetMax);

% plot results
t = [0:length(uvalid)-1]';
f3=figure('Name', 'Exact Simulation');
plotgp(f3,t,yvalid, y, sqrt(s2));