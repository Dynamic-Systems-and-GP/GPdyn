%% detect_fault_validity
function [index_validity, index_validity_window] = detect_fault_validity...
    (X, input, target, NW, var_residuals)

%% Syntax
%  function [index_validity, index_validity_window] = detect_fault_validity...
%    (X, input, target, NW, var_residuals)


%% Description
% Function, which estimates the degree at which we can trust the results of 
% the function DETECT_FAULT. To measure the quality of the identified model
% particular region, the GP model prediction variances are used.  
% Note: 
% (1) Function is not compatible with the current version of the GPML toolbox. 
% (2) The simulation methods which propagate the variance are recommended, 
% eg. simul00mcmc etc. 
% [1] Dj.Juricic, J.Kocijan, Fault detection based on Gaussian process model,
% I.Troch, F.Breitenecker (eds.), Proceedings of the 5th Vienna
% Symposium on Mathematical Modelling – MathMod, Wien, 2006
% 
% Inputs: 
% X, input, target .. standard GP model inputs 
% NW .. window size 
% var_residuals .. vector of predicted variances
% var_residuals .. [v(1) .. v(n)], where w(k) = var(eps(k))  predicted
%       variance associated with yhat(k), i.e. GP prediction at time step k
% NW .. size of window, used for calculating the index and deciding about fault
% c_max .. maximal significance level still reflecting in deciding for H0
%   (no fault)
% Outputs: 
% index_validity: validity index for time steps Ik 
%   (Ik=0-> not valid, Ik=1 -> valid)
% index_Wvalidity: validity index Ik averaged over window NW 

%% See Also
% DETECT_FAULT, SIMUL02MCMC


% set to 1 while debugging and testing
flag_plot = 0;


KK = []; 

delta_max = 0; 
% we determin delta_max 
for kk=1:size(input,1)
    test = input(kk,:); 
    [dumm,kandidat] = gpSD00(X, input, target, targetvariance, derivinput, derivtarget, derivvariance, test);
    KK = [KK; kandidat]; 
    if(kandidat > delta_max)
        delta_max = kandidat; 
    end 
end 


figure(45)
plot(KK); 
hold on 
plot(delta_max*ones(size(kandidat))); 
hold off 


% validity index 
Ik = []; 


for kk=1:length(s2)
    
    % option 1 (not used): determine the difference between the "mean
    % input" to the GP model and the training input (neede additional input
    % - input regressors to GP model for all time steps) 
    if(0)  
        test = regressors(kk,:); 
        [dumm,delta] =  gpSD00(X, input, target, targetvariance,...
        derivinput, derivtarget, derivvariance, test);
        %delta = delta + exp(2*X(end)); 
    % option 2: distance between test and training inputs is estimated from
    % the predicted variances 
    else
        delta = s2(kk)-exp(2*X(end)); 
    end 

   
    if(delta<=delta_max)
        Ik(kk) = 1; 
    else
        Ik(kk) = delta_max/delta; 
    end 
end 

for kk=1:length(Ik)
    if(kk<NW)
        IkW(kk) = mean(Ik(1:kk)); 
    else
        IkW(kk) = mean(Ik(kk-NW+1:kk)); 
    end 
end 

index_validity = Ik'; 
index_validity_window = IkW'; 



return;

