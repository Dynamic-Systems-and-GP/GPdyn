%% mcmc_test_pdfs
function [] = mcmc_test_pdfs(MM,VV,desired_steps)

%% Syntax
% function [] = mcmc_test_pdfs(MM,VV,desired_steps)

%% Description
% Function test and plots the Gaussian mixture given my MM and VV in
% desired simulation steps 
% Inputs: 
% MM .. mean values of predicted distributions, matrix ksteps times Nsamples
% VV .. associated variances of predicted distributions 
% desired_steps .. vector with the indices of steps, where we want to test
% the distributions 
% Outputs: 
% .. plots the figures of the distributions in the desired steps 

%% See Also
%
% SIMUL00EXACT



% XX,YY ... matrices for storing the GP model's pdfs
% XX(:,k) ... yvalues at time step k
% YY(:,k) ... density values at time step k and coresponding y values in XX(:,k)
% sum(YY(:,k)) = 1 for all k-s

% variables 
MM = MM'; 
VV = VV'; 

Nsamples = size(MM,1); 
mu = mean(MM);
s2 = mean(VV) + mean((MM-repmat(mu,Nsamples,1)).^2);
mu = mu';
s2 = s2';

% "calculation"
STD = sqrt(VV);
Xwide = 200; % resolution = no. of samples between xmin and xmax
Xmin = min(MM-4*STD);
Xmax = max(MM+4*STD);
dX = (Xmax-Xmin)/(Xwide-1);

XX = zeros(Xwide,size(VV,2));
YY = zeros(size(XX));


% for desired steps
for jj = 1:length(desired_steps)
    kk = desired_steps(jj);
    XX(:,kk) = Xmin(kk):dX(kk):Xmax(kk);
    % ... and for every realisation
    ytemp = zeros(size(YY,1),1);

    for ii=1:Nsamples
        ytemp = normpdf(XX(:,kk), MM(ii,kk), STD(ii,kk));
        YY(:,kk) = YY(:,kk) + ytemp;
    end
    YY(:,kk) = YY(:,kk)/Nsamples;
    
    yapprox = normpdf(XX(:,kk),mu(kk),sqrt(s2(kk)));

    % plot pdfs 
    figure(10000+kk);

    plot(XX(:,kk),YY(:,kk),'Color',[0.6 0.6 0.6],'LineWidth',4);
    hold on;
    plot(XX(:,kk),yapprox,'LineStyle','--','Color',[0 0 1],'LineWidth',4);
     for ii=1:Nsamples
        ytemp = normpdf(XX(:,kk), MM(ii,kk), STD(ii,kk));
        plot(XX(:,kk),ytemp,'Color',[0.9 0.9 0.9]);
    end    
    plot(XX(:,kk),YY(:,kk),'Color',[0.6  0.6 0.6],'LineWidth',4);
    plot(XX(:,kk),yapprox,'LineStyle','--','Color',[0 0 1],'LineWidth',4);
    hold off;
    grid; 
    legend('true','GP approx','samples')


    title(strcat(['test mcmc pdfs, simul step=', num2str(kk)]));


    disp('press sth to continue'); 
    pause;

    
end
