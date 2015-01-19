%% mcmc_getsamplesgaussmix
function [x,y,indexes,samples,xsamples,fwarning] = mcmc_getsamplesgaussmix(x,Mi,Sig,N)

%% Syntax
% function [x,y,indexes,samples,xsamples,fwarning] = mcmc_getsamplesgaussmix(x,Mi,Sig,N)

%% Description
% Function returns the samples from the mixtures of Gaussian distributions,
% given by Mi and Sig and also estimates the quality of the sampled
% distribution. First the cumulative distribution is estimated, then the
% inverse transform method with the use of Matlab random ganerator is used
% to obtain individual samples. 
% Inputs: 
% x .. partitioned input space
% Mi .. vector of Gaussian mixture means
% Sig .. vectors of Gaussian mixture standard deviations 
% N ... number of returned samples
% Outputs: 
% x .. partitioned inpout space 
% y .. probabitily density function (pdf) of Gaussian mixture 
% indexes .. indexes of chosen samples, should be applied to x to get the
%   values 
% samples .. number of samples at particular elements of x - for testing
%   the quality of the approximated distribution
% xsamples .. indexes applied to x (samples) 
% fwarning .. warning flag, risen if less than 98% of distribution is
%   covered with samples 
% 
% Function uses rand (internal Matlab command) and mcmc_getgaussmix,
% defined at the end of the function. It is used by the function
% simul00mcmcm. 

%% See Also
%
% RAND, SIMUL00MCMC


[x,y,cumy] = mcmc_getgaussmix(x, Mi, Sig);
% to get "real" sum(y) (==1) one must divide y with (1/dx) i.e.
% yreal = y / (1/dx) = y*dx i.e.
% y = y*(x(2)-x(1));

fwarning = 0; 
% if too big part of distribution is left out we print the warning
RATIO = 0.98;
ratio = (x(2)-x(1))*sum(y);
if(ratio<RATIO)
    warning(['at least ',num2str(100*RATIO), '% of the complete distribution should be sampled, only ', num2str(100*ratio,3), '% is.']);
    fwarning = 1; 
    %keyboard; 
end


%tic
samples = zeros(size(x));
xi = rand(N,1);
for ii=1:N
    % inverse transform method - check where uniform random number "falls"
    % into cumulative distribution

    jj = find(xi(ii)<cumy,1);

    % upading vector of indexes - for sampling
    indexes(ii) = jj;
    % increasing the number of samples at corresponding x by one - for
    % analysis (to check whether samples are describing desired
    % distribution)
    samples(jj) = samples(jj)+1;
end
samples = samples * sum(y)/N;
%toc


% if needed samples are already applied to x and the outputs are ready for
% use
x = [x(1)-x(2)+x(1); x];  % CC1 if min(indexes)==x(1)
if (nargout>=5)
    xsamples = (x(indexes)+x(indexes+1))/2;
end
x = x(2:end); % check comment CC1 4 rows up

if (fwarning==1)
    dx = x(2)-x(1);
    figure(11);
    plot(x,samples*dx);
    hold on;
    plot(x,y*dx,'LineWidth',1,'Color',[0.7 0.7 0.7]);
    hold off;
    legend('samples','real mixture');
    %keyboard; 
end


return;




% ***********************************************************************
% ***********************************************************************

function [x,y,cumy] = mcmc_getgaussmix(x, Mi, Sig);
% Subfunction, which returns the probability and cumulative distribution 
% for the (with Mi and Sig) defined Gaussian mixture.   


y = zeros(size(x));
for i=1:length(Mi)
    yi = normpdf(x, Mi(i), Sig(i));
    y = y + yi; 
end
y = y/length(Mi);

if(nargout==3)
    cumy = cumsum(y);
    cumy = cumy/cumy(end);
end

