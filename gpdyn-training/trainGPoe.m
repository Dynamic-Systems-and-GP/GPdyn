function [hyp, flogtheta, i,simy,simse2] = trainGPoe(hyp, inf, mean, cov, lik, input, target, lag, simf, Nsamples)
% Function for the optimisation (training) of GP model hyperparameters
%
%% Syntax
% [hyp, flogtheta, i] = trainGPoe(hyp, inf, mean, cov, lik, input, target,
%                                 simf, lag, Nsamples);
%
%% Description
% Function for the optimisation (training) of GP model hyperparameters
% based on the training data via maximum marginal likelihood. 
% Uses routines gp and minimize.
% Based on the work of C.E.Rasmussen. 
% 
% Input:
% * hyp      ... the structure of initial hyperparameters
% * inf      ... the function specifying the inference method 
% * cov      ... the prior covariance function (see below)
% * mean     ... the prior mean function
% * lik      ... the likelihood function
% * input    ... the input part of the training data,  NxD matrix
% * target   ... the output part of the training data (ie. target), Nx1 vector
% * simf     ... the function handle of simulation function (e.g. @simulGPmc)
% * lag      ... the order of the model (number of used lagged outputs)
% * Nsamples ... the number of samples for MCMC simulation (optional)
% * minf     ... the function handle of the minimization method to be used 
%                (optional, default=@minimize_and_return). The minimization 
%				  method should return the following values:
%				  [X, fX, i,fout] = minimize(X, f, length, P1, P2, P3, ... ),
%				  where fout is the additional argument (along fitness and gradient)
%				  of the evaulation function f
%
% Output: 
% * hyp       ... optimized hyperparameters 
% * flogtheta ... the minus log likelihood for the different runs (init. to 0)
% * i         ... the number of iterations needed for the last optimization
%
% Examples:
% demo_example_gp_training.m
%
% See Also:
% gp, minimize, covFunctions, trainlgmp
%

if(nargin < 11)
  minf = @minimize;
end
if (nargin < 10)
  Nsamples = 100; % default value of Nsamples for MC simulation
end
if (nargin < 9)
  simf = @simulGPnaive; %choose the default simulation routine
end

if (nargin < 8)
  warning(['The number of dynamical system model order was not provided. ' ...
  'Guessing from the input data...']);
  lag=size(input,2)/2;
  assert(mod(lag,1)==0);
  fprintf('Guess: lag=%d\n',lag);
end

if (nargin < 7)
  error('Too few parameters are given.');
end


MIN_DIFF = 0.002; 
flogtheta = [0,-1/eps]; 

% we reconstruct the input signal from the regressor matrix as
%  the two following commented lines show:
%%N = length(input,1) +    lag-1
%%[u(1) ... u(N-lag)]';	[u(N-lag+1) ...  u(N)]
 
u=[input(:,lag+1)      ;	input(end,lag+2:end);NaN];

simulated_y=[input(1:lag,1);target];

while (abs(flogtheta(end) - flogtheta(end-1))>MIN_DIFF)         
	disp(' '); 
	disp(strcat(['delta flogtheta: ', num2str(abs(flogtheta(end) - flogtheta(end-1)))])); 
	disp(' ')
	[hyp, flogthetatmp,i] = feval(minf, hyp, @simLL, -10,...
				inf, mean, cov, lik, u, simulated_y, target, lag, simf);
	
	%in the last iteration of minimisation algorithm the simulation is evaluated again for obtaining
	% the simulated response 'simulated_y'
	disp('saving simulated response.')
    [~,~,simulated_y] = simLL(hyp,inf, mean, cov, lik, u, simulated_y, target, lag, simf);
    
    if isempty(flogthetatmp) % no improvement: at minimum
        disp('oops')
        break
    end
	flogtheta = [flogtheta flogthetatmp(end)];
end


function [nlZ,dnlZ,simulated_y]=simLL(hyp, inf, mean, cov, lik, u, simulated_y, target,lag,simf)
	% simLL simulates the response and forms the regressor vectors from the response. Despite the response signal in regressor vectors is the simulated one, the regressands are always the original measurements.
	%
	% log-likelihood is calculated from the gp model with input data from simulation (in regressor matrix) and measurements (as regressands).
	%
	% simulated_y is a NxS or 2NxS matrix where N is the signal length and S is the number of past iterations. The values from N+1 to 2N represent the predicted variances.
	%
	N=size(target,1)+lag;
	if size(simulated_y,1)==N
		simulated_y=[simulated_y;simulated_y*NaN];
	end
	
	simulated_input0 = construct(lag, u, simulated_y(1:N,end)); 
	
	[y_mu, se2] = feval(simf,hyp, inf, mean, cov, lik, simulated_input0, ...
	target, simulated_input0, lag);
	
	%y_mu_full combines measured initial values and simulated response:
	simulated_y(1:N,end+1)=[simulated_y(1:lag,end);y_mu(:)];
	simulated_y(N+1:2*N,end)=[zeros(lag,1);se2(:)];
	
	simulated_input1 = construct(lag, u, simulated_y(1:N,end)); 
	[nlZ,dnlZ]=gp(hyp,inf,mean,cov,lik,simulated_input1,target);
end

simy=simulated_y(1:end/2,:);
simse2=simulated_y(end/2:end,:);
end


