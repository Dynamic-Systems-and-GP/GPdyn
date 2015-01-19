function [nInput,inputMin,inputMax, nTarget, targetMin, targetMax] = preNorm(input, target, inputMax)
%Preprocesses data so that minimum is -1 and maximum is 1.
%  
%% Syntax
%   [nInput,inputMin,inputMax] = preNorm(input);
%   [nInput,inputMin,inputMax,nTarget,targetMin,targetMax] = preNorm(input, target);
%   [nInput] = preNorm(input,inputMin,inputMax);
% 
% 
%% Description
% Function normalizes the inputs so that they fall in the interval [-1,1].
% Algorithm:  nInput = 2*(input-inputMin)/(inputMax-inputMin) - 1; If tde
% data needs to be normalzed using previously coputed minimums and
% maximums, we use: nInput = preNorm(input,inputMin,inputMax);
% 
% Input:
% * input ... n x D input matrix
% * target/inputMin ... target vector/if there are 3 input arguments it
%                       will be considered as a vector of minimums
% * inputMax ... maximums
% 
% Output:
% * nInput ... n x D normalized input matrix
% * inputMin ... row vector containing minimums for each dimension
% * inputMax ... row vector containing maximums for each dimension
% * nTarget ... normalized target vector
% * targetMin ... target minimum
% * targetMax ... target maximum
% 
% See also:
% postNorm, postNormVar

%%
% * Written by Tomaž Šuštar, January 2012



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
  warning('Some maximums and minimums are equal. Those inputs won''t be transformed.');
  inputMin0 = inputMin.*notequal - 1*isequal; % where equal set minimums to -1
  inputMax0 = inputMax.*notequal + 1*isequal; % and maximums to +1 so the data will not be transformed
else
  inputMin0 = inputMin;
  inputMax0 = inputMax;
end

nInput = 2*(input-repmat(inputMin0,n,1))./repmat((inputMax0-inputMin0),n,1) - 1; % normalize



