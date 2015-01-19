function [data, target] = postNorm(nInput,inputMin,inputMax, nTarget, targetMin, targetMax)
% Postprocesses data which has been preprocessed by preNorm.
%  
%% Syntax
%    [data] = postNorm(nInput,inputMin,inputMax);
%    [data,target] = postNorm(nInput,inputMin,inputMax,nTarget,targetMin,targetMax);
%
%% Description  
% This function postprocesses the training set which was preprocessed by
% preNorm function. It converts the data back into unnormalized units.
% Algorithm: data = 0.5(nInput+1)*(inputMax-inputMin) + inputMin;
%  
% Input:
% * nInput ... n x D normalized input matrix
% * inputMin ... row vector containing minimums for each dimension
% * inputMax ... row vector containing maximums for each dimension
% * nTarget ... normalized target vector
% * targetMin ... target minimum
% * targetMax ... target maximum
% 
% Output:
% * data ... unnormalized data
% * target ... target vector  
%
% See also:
% preNorm, postNormVar
%%
% * Written by Tomaž Šuštar, January 2012

data = unnormalize(nInput,inputMin,inputMax);

if nargin==6
  target = unnormalize(nTarget, targetMin, targetMax);
end

function output = unnormalize(nInput,inputMin,inputMax)

[n, D] = size(nInput); % n - number of mesurements, D - input space dimenson

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

output = (nInput+1)/2 .* repmat((inputMax0-inputMin0),n,1) + repmat(inputMin0,n,1); % unnormalize



