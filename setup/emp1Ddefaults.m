% default inputs for emp 2D jobs

inputs.exefile = 'emp1d';
inputs.exedir = '/shared/users/ram80/empcodes/emp1/';

inputs.submitjob = 1;  % set to zero just to test setup

% do you want to include the ionosphere?
inputs.doionosphere = 1;
% do you want to calculate ionosphere changes (Ne, etc)?
inputs.doioniz = 1;
% number of times to write to output arrays - evenly distributed
inputs.numfiles = 30;
% highest altitude to consider. Notice everything in meters!
inputs.maxalt = 300e3;

% initial dx from ground up
inputs.dx1 = 200;
% dx at higher altitudes
inputs.dx2 = 200;
% altitude at which to change to smaller dx. If they are the same, it
% doesn't matter.
inputs.stepalt = 70e3;

% time step: use 1e-7 for D-region (good up to 150 km, grids >= 100 m)
inputs.dt = 1e-7;

% lightning inputs
inputs.E0 = 100; % source amplitude in V/m
inputs.taur = 10e-6;
inputs.tauf = 30e-6;

% filtering for the lightning input
inputs.fcut = 300e3;

% magnetic field
inputs.Bmag = [-35000 0 35000] * 1e-9;

% cluster and number of nodes
inputs.cluster = 'batchnew';
inputs.numnodes = '8';

% probe altitudes
inputs.probealt = [100 150 200 250 301 350] * 1e3;
