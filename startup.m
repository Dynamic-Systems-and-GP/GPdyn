% startup script to make Matlab aware of the GpDyn package

disp(['executing gpdyn startup script...']);

OCT = exist('OCTAVE_VERSION') ~= 0;           % check if we run Matlab or Octave

me = mfilename;                                            % what is my filename
mydir_gpdyn = which(me); mydir_gpdyn = mydir_gpdyn(1:end-2-numel(me));        % where am I located
if OCT && numel(mydir_gpdyn)==2
  if strcmp(mydir_gpdyn,'./'), mydir_gpdyn = [pwd,mydir_gpdyn(2:end)]; end
end                 % OCTAVE 3.0.x relative, MATLAB and newer have absolute path

addpath(mydir_gpdyn(1:end-1))                                  % add dirs to path
addpath([mydir_gpdyn,'gpdyn-gp-evaluation'])
addpath([mydir_gpdyn,'gpdyn-lmgp-evaluation'])
addpath([mydir_gpdyn,'gpdyn-training'])
addpath([mydir_gpdyn,'gpdyn-utilities'])
addpath([mydir_gpdyn,'gpml-matlab'])

run([mydir_gpdyn,'gpml-matlab/startup.m'])          % execute gpml startup script

clear me mydir OCT
