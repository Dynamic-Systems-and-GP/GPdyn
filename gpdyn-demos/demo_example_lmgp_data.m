% demo_example_lmgp_data - demonstration script file for LMGP model identification.  
%
%% Description
% Demo to present how to obtain and compose data for the identification
% of the GP model with incorporated local models (LMGP model). 
% 
% There are two options how to obtain local model data: 
% a. analyticaly (set fIdentifyLm to 0) 
% b. with identification (set fIdentifyLm to 1) 
% 
% It can be used only with Gaussian covariance function and
% with white noise model (sum of covSEard and covNoise). 
%
% See Also
% example.m, demo_example_LMGP_training.m, demo_example_LMGP_simulation.m,
% example_LM_ident.m, example_derivative.m
%
%% 
% Written by K. Azman, 2007


clear;
global flag_LM_data_ident
% LM generation data 
flag_LM_data_ident = 0; % 0=calculated, 1=identified using IV

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
    
    filename = 'example_data_LMGP_ident';
    
else
    % *** local model's parameters calculated ***
    [eq.dfdy, eq.dfdu, eq.Y] = demo_example_derivative(eq.U); 
    
    filename = 'example_data_LMGP_anal';
   
end

save(filename,'eq')
    
return 