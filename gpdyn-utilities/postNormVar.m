function s2 = postNormVar(s2n,targetMin,targetMax) 
% Postprocesses predicted variance for data which has been preprocessed by 
% preNorm.
%
%% Syntax
%  s2 = postNormVar(s2n,min,max) 
%
% Input: 
% * s2n       ... the normalized predicted variance 
% * targetMin ... the target minimum
% * targetMax ... the target maximum
%
% Output: 
% * s2 ... the postprocessed predicted variance 
%
%% See Also
% preNorm, postNorm 

stdn = sqrt(s2n); 
std = stdn*(targetMin-targetMax)/2; 
s2 = std.^2; 

return 