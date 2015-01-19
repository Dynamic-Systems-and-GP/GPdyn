% startup script to make Matlab aware of the GpDyn package

disp(['executing gpdyn startup script...']);

me = mfilename;                                       % what is my filename
mydir = which(me); mydir = mydir(1:end-2-numel(me));   % where am I located

addpath(mydir(1:end-1))                                  % add dirs to path
addpath([mydir,'gpdyn-gp-evaluation'])
addpath([mydir,'gpdyn-lmgp-evaluation'])
addpath([mydir,'gpdyn-training'])
addpath([mydir,'gpdyn-utilities'])
addpath([mydir,'gpml-matlab'])

run([mydir,'gpml-matlab/startup.m'])          % execute gpml startup script

clear me mydir OCT
