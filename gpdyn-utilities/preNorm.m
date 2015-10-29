function [nInput,inputMin,inputMax, nTarget, targetMin, targetMax] = preNorm(input, target, inputMax)
% Preprocesses data so that minimum is -1 and maximum is 1.
%  
%% Syntax
%   [nInput,inputMin,inputMax] = preNorm(input);
%   [nInput,inputMin,inputMax,nTarget,targetMin,targetMax] = preNorm(input, target);
%   [nInput] = preNorm(input,inputMin,inputMax);
% 
% 
%% Description
% Function normalizes inputs so that they fall in the interval [-1,1].
% Algorithm:  nInput = 2*(input-inputMin)/(inputMax-inputMin) - 1; If the
% data needs to be normalzed using previously computed minimums and
% maximums, we use: nInput = preNorm(input,inputMin,inputMax);
% 
% Input:
% * input           ... the n x D input matrix
% * target/inputMin ... the target vector/if there are 3 input arguments it
%                       will be considered as a vector of minimums
% * inputMax        ... maximums
% 
% Output:
% * nInput      ... the n x D normalized input matrix
% * inputMin    ... the row vector containing minimums for each dimension
% * inputMax    ... the row vector containing maximums for each dimension
% * nTarget     ... the normalized target vector
% * targetMin   ... the target minimum
% * targetMax   ... the target maximum
% 
% See also:
% postNorm, postNormVar

%% Signature
% * Written by Tomaz Sustar, January 2012



if nargin > 3
  error('Wrong number of arguments.');
end

if nargin==3 % normalize targets if given
  inputMin = target;
  [nInput,inputMin,inputMax] = normalize(input,inputMin,inputMax); %normalize input
else
  [nInput,inputMin,inputMax] = normalize(input); %normalize input
end


if nargin==2 % normalize targets if given
  [nTarget, targetMin, targetMax] = normalize(target); 
end



function [nInput,inputMin,inputMax] = normalize(input, inputMin, inputMax)
[n, D] = size(input); % n - number of mesurements, D - input space dimenson 

if nargin ~= 3
  inputMin = min(input); % row vector of minimiums 
  inputMax = max(input); % row vector of maximums
end

isequal = inputMin==inputMax;
notequal = ~isequal;
if sum(isequal) ~= 0
  warning('Some maximums and minimums are equal. Those inputs will not be transformed.');
  inputMin0 = inputMin.*notequal - 1*isequal; % where equal set minimums to -1
  inputMax0 = inputMax.*notequal + 1*isequal; % and maximums to +1 so the data will not be transformed
else
  inputMin0 = inputMin;
  inputMax0 = inputMax;
end

nInput = 2*(input-repmat(inputMin0,n,1))./repmat((inputMax0-inputMin0),n,1) - 1; % normalize



