%% postmnmxvar
function s2 = postmnmxvar(s2n,min,max) 

%% Syntax
%  function s2 = postmnmxvar(s2n,min,max) 

%% Description
% Postprocesses predicted variance for data which has been preprocessed by 
% PREMNMX.
% Inputs: 
% s2n .. normalized predicted variance 
% min .. minimum of preprocessed data 
% max .. maximum of preprocessed data 
% Output: 
% s2 .. postprocessed predicted variance 

%% See Also
% POSTMNMX, PREMNMMX 

stdn = sqrt(s2n); 
std = stdn*(max-min)/2; 
s2 = std.^2; 

return 