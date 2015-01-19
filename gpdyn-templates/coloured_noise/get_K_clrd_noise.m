%% get_K_clrd_noise
function [Cn, Cn_vec] = get_K_clrd_noise(NN, a, b, noisevar)

%% Syntax
%  function [Cn, Cn_vec] = get_K_clrd_noise(NN, a, b, noisevar)

%% Description
% Function returns the noise part of covariance matrix, when the 
% coloured noise on the output of the system is known. 
% This covariance matrix should be added to its function part. The
% function folows the derivation in [1] and differs a bit from [2]. 
% Idea:  
% AR part (ai-s) : 
% first na elements of the covariance matrix calculated using Yule-Walker
% equations, the rest of them (na+1 .. NN) computed iteratively 
% MA part (bi-s) : 
% the elements of covariance matrix up to bn are calculated, the rest of
% them equals 0 
% Notes:  
% (1) K(tau) = K(i,j), where |i-j|=|j-i|=tau 
% (2) Function is not compatible with the current version of the GPML 
% toolbox. 
% (3) Function uses the function TOEPLITZ. 
% 
% Inputs: 
% NN .. size of the covariance matrix (number of training data) 
% a ... vector of AR elements of correlated-noise (1 a_1 a_2 .. a_na)
% b ... vector of MA elements of correlated-noise (b_0 b_1 b_2 .. b_nb)
%    usually b_0=0
% noisevar ... noise variance sigma_e^2 
% Outputs: 
% Cn .. noise part of the covariance matrix 
% Cn_vec .. elements of the covariance matrix in a vector 

%% See Also
% [1] K. Azman, J. Kocijan, Identifikacija dinamiènega sistema z znanim modelom
% šuma z modelom na osnovi Gaussovih procesov, n: B. Zajc, A. Trost (eds.)
% Zbornik petnajste mednarodne Elektrotehniške in raèunalniške konference 
% ERK 2006, 25.-27.09., Portorož, Slovenija, 2006, Vol. A, pp. 289--292, 
% 2006.  (in Slovene) 
% [2] R. Murray-Smith and A. Girard. Gaussian Process Priors with ARMA 
% noise models. In Proceedings: Irish Signals and Systems Conference, 
% 147-153, Maynooth, 2001.


flag_predznak = -1; % -1: Azman 

% cut first elements (a0=1 and b0=0)
%a = a(1:end);   %  a was equal a=[ 1 a1 a2 .. a_Nb]
%b = b(1:end);   %  b was equal b=[b0 b1 b2 .. b_Nb]
a0 = a(1); 
a = a(2:end); 
b0 = b(1);
b = b(2:end);
Na = length(a);
Nb = length(b); 

% horizontal llines 
if(size(a,1)>size(a,2))
    a=a';
end
if(size(b,1)>size(b,2))
    b=b';
end

if(Na > Nb)
    warning('AR and MA model have different lengths. zeros added to b'); 
    b(Nb+1:Na) = zeros(1, Na-Nb); 
elseif(Nb>Na) 
    warning('AR and MA model have different lengths. zeros added to a'); 
    a(Na+1:Nb) = zeros(1, Nb-Na); 
end 



% generally the size of covariance matrix > dynamics order 
if(NN <= max(Na,Nb))
    error('order of noise model to high or the size of training matrix not big enough');
end


% ************************************************************************
% MA part 
% no MA part 
if(Nb==0)
    MA0 = noisevar*b0^2;
    MA = zeros(1,NN-1);
else
    % diagonal elements 
    MA0 = noisevar*(b0^2+sum(b.*b));
    % rest 
    for ii=1:Nb
        vec1 = [b(ii:end)];
        vec2 = [b0 b(1:end-ii)];        
        MA(ii)= sum(vec1.*vec2);
    end 
    MA = noisevar*MA; 
    %  MA part of the noise model represented with MA0 and vector MA
end 



% ************************************************************************
% AR part 
if(Na>0)

    a = [a0 -flag_predznak*a]; 
% index matrix  
Ain = toeplitz(1:Na+1); 
% generation of YW matrix 
A = zeros(Na+1); 
for ii=1:Na+1
    for jj=1:Na+1
        stolpec = Ain(ii,jj); 
        A(ii,stolpec) = A(ii,stolpec)+ a(jj); 
    end 
end 

% A*AR = B, where: 
% A upper matrix 
% AR covariance vector (we have to calculate it): C_0..C_Na
% B vector of MA part 
AR = inv(A)*[MA0 MA]'; 
AR0 = AR(1); 
AR = [AR(2:end)]'; 
a = a(2:end); 

AR0 = MA0; 
AR = MA; 

% the rest of the covariance matrix elements 
if (Na>0)
    for ii=Na+1:NN-1 % one less than in NN, as one already in MA0
        AR(ii)= flag_predznak*sum(a.*AR(ii-1:-1:ii-Na));
    end
end 



% ************************************************************************
% combine the values of both vectors into covariance matrix Cnoise
Cn_vec = [AR0 AR];
Cn_vec = [Cn_vec zeros(1,NN-length(Cn_vec))]; 

% building the matrix from vector 
Cn = toeplitz(Cn_vec);









