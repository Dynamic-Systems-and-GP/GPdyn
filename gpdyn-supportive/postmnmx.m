function [p,t] = postmnmx(pn,minp,maxp,tn,mint,maxt)
%POSTMNMX Postprocesses data which has been preprocessed by PREMNMX.
%  
%  Syntax
%
%    [p,t] = postmnmx(pn,minp,maxp,tn,mint,maxt)
%     [p] = postmnmx(pn,minp,maxp)
%
%  Description
%  
%    POSTMNMX postprocesses the network training
%    set which was preprocessed by PREMNMX.  It converts
%     the data back into unnormalized units.
%  
%    POSTMNMX takes these inputs,
%       PN  - RxQ matrix of normalized input vectors
%       MINP- Rx1 vector containing minimums for each P
%       MAXP- Rx1 vector containing maximums for each P
%       TN  - SxQ matrix of normalized target vectors
%       MINT- Sx1 vector containing minimums for each T
%       MAXT- Sx1 vector containing maximums for each T
%    and returns,
%      P - RxQ matrix of input (column) vectors.
%       T - SxQ matrix of target vectors.
%    
%  Examples
%
%    In this example we normalize a set of training data with
%     PREMNMX, create and train a network using the normalized
%     data, simulate the network, unnormalize the output of the
%     network using POSTMNMX, and perform a linear regression between 
%     the network outputs (unnormalized) and the targets to check the
%     quality of the network training.
%  
%      p = [-0.92 0.73 -0.47 0.74 0.29; -0.08 0.86 -0.67 -0.52 0.93];
%       t = [-0.08 3.4 -0.82 0.69 3.1];
%      [pn,minp,maxp,tn,mint,maxt] = premnmx(p,t);
%       net = newff(minmax(pn),[5 1],{'tansig' 'purelin'},'trainlm');
%       net = train(net,pn,tn);
%       an = sim(net,pn);
%       [a] = postmnmx(an,mint,maxt);
%       [m,b,r] = postreg(a,t);
%
%  Algorithm
%
%     p = 0.5(pn+1)*(maxp-minp) + minp;
%
%  See also PREMNMX, PREPCA, POSTSTD.

% Copyright 1992-2002 The MathWorks, Inc.
% $Revision: 1.7 $

[R,Q]=size(pn);
oneQ = ones(1,Q);

equal = minp==maxp;
nequal = ~equal;
if sum(equal) ~= 0
  warning('Some maximums and minimums are equal. Those inputs won''t be transformed.');
  minp0 = minp.*nequal - 1*equal;
  maxp0 = maxp.*nequal + 1*equal;
else
  minp0 = minp;
  maxp0 = maxp;
end

p = (pn+1)/2;
p = p.*((maxp0-minp0)*oneQ) + minp0*oneQ;

if nargin==6
  equal = mint==maxt;
  nequal = ~equal;
  if sum(equal) ~= 0
    warning('Some maximums and minimums are equal. Those targets won''t be transformed.');
    mint0 = mint.*nequal - 1*equal;
    maxt0 = maxt.*nequal + 1*equal;
  else
    mint0 = mint;
    maxt0 = maxt;
  end
  t = (tn+1)/2;
  t = t.*((maxt0-mint0)*oneQ) + mint0*oneQ;
end
