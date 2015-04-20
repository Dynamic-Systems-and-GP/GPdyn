function [] = mc_test_pdfs(MM,VV,desired_steps)
% Function tests and plots the Gaussian mixture in the desired simulation
% steps. 
%
%% Syntax
% mc_test_pdfs(MM,VV,desired_steps)
%
%% Description
% Function tests and plots the Gaussian mixture described my vector of Gaussian
% means and variances in the  desired simulation steps.
%
% Input: 
% * MM   ... matrix of mean values of Gaussian components, the matrix of dimensions 
%            ksteps x Nsamples
% * VV   ... matrix of associated variances of Gaussian components 
% * desired_steps ... vector with indices of steps, where we want to test
%            the distributions 
%
% See Also:
% demo_example_gp_simulation
% 
% Examples: 
% demo_example_gp_simulation
% 

%% 
% * Written by K. Azman, 2007.
%

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

% XX,YY ... matrices for storing the GP model's pdfs
% XX(:,k) ... yvalues at time step k
% YY(:,k) ... density values at time step k and coresponding y values in XX(:,k)
% sum(YY(:,k)) = 1 for all k-s
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


    title(strcat(['test MC pdfs, simul step=', num2str(kk)]));


    disp('press sth to continue'); 
    pause;

    
end
