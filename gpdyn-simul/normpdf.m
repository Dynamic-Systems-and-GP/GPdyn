function p = normpdf ( x, mean, sigma )
% NORMPDF.M calculates and plots a
% normal probability density function.

%===============================================
% mean   = mean of the data vector, x
% sigma  = standard deviation of data vector, x
% x      = data vector over which to find p(x)
% p(x)   = probability density function
%===============================================

normalizer = 1/(sigma*sqrt(2*pi));
p = normalizer*exp( -0.5*( (x - mean)/sigma ).^2 );

