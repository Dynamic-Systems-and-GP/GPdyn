function s2 = postNormVar(s2n,targetMin,targetMax) 
% Postprocesses predicted variance for data which has been preprocessed by 
% preNorm.
%
%% Syntax
%  s2 = postNormVar(s2n,min,max) 
%
% Input: 
% * s2n ... normalized predicted variance 
% * targetMin ... target minimum
% * targetMax ... target maximum
%
% Output: 
% * s2 ... postprocessed predicted variance 
%
%% See Also
% preNorm, postNorm 

stdn = sqrt(s2n); 
std = stdn*(targetMin-targetMax)/2; 
s2 = std.^2; 

return 