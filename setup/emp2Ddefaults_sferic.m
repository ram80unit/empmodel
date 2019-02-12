% default inputs for emp 2D jobs

inputs.exefile = 'emp2d';
inputs.exedir = '/shared/users/ram80/empcodes/emp2/';

inputs.submitjob = 1;  % set to zero just to test setup
inputs.savefields = 1; % set to zero if you don't need the large output fields

% what planet?
inputs.Re = 6370000;

% do you want to use the PML boundary?
inputs.dopml_top = 0;
inputs.dopml_wall = 1;
% do you want to include the ionosphere?
inputs.doionosphere = 0;
% do you want to calculate ionosphere changes (Ne, etc)?
inputs.doioniz = 0;
% do you want to integrate and spit out the elve movie?
inputs.doelve = 0;
% do you want to do associative detachment?
inputs.dodetach = 0;
% number of times to write to output arrays - evenly distributed
inputs.numfiles = 30;
% highest altitude to consider. Notice everything in meters!
inputs.maxalt = 85e3;
% perfectly conducting ground, or real ground? 
inputs.pecground = 0;

% dr below surface (ground)
inputs.dr0 = 20;
%number of ground cells
inputs.nground = 10;
% initial dr from ground up
inputs.dr1 = 100;
% dr at higher altitudes
inputs.dr2 = 100;
% delta in range direction (used to calculate dtheta)
inputs.drange = 100;
% altitude at which to change to smaller dr. If they are the same, it
% doesn't matter.
inputs.stepalt = 70e3;

% lightning location, azimuth of interest, and range
inputs.Trlat = 34.25;
inputs.Trlon = -88.97;
inputs.range = 510e3;
inputs.az = 298;

% time step: use 1e-7 for D-region (good up to 150 km, grids >= 100 m)
inputs.dt = 1e-7;

% conductivities: make them nonzero if you want
inputs.sig = 0;
inputs.sigm = 0;

% camera location, distance from source along ground and altitude
inputs.camdist = 400e3;
inputs.camalt = 0;
inputs.camelev = 9;
inputs.camfov = [36 18];  %% left-right, up-down
inputs.numpixels = [128 64];   %% camera pixels
inputs.cameratype = 'camera';
inputs.elvesteps = 1000;
% transmitter? venus?
inputs.dotransmitter = 0;
inputs.dovenus = 0;
inputs.txf0 = 20e3;

%decimate outputs before writing. Useful for 100m resolution. Default = 1
inputs.decfactor = 5;

% set up probe points. Define locations in km, then determine grid values
inputs.probealt = [0.5 0.5 0.5 0.5 0.5] * 1e3;
inputs.proberange = [100 200 300 400 500] * 1e3;

% lightning inputs
inputs.I0 = 100e3; %e3; % peak current in amperes:
inputs.Ic = 0; %2e3;   % continuing current!
inputs.sourcealt = 10e3;
inputs.taur = 7e-6;
inputs.tauf = 50e-6;
inputs.rsspeed = 2.5e8;
inputs.decaytype = 1;
% 0 is TL, 1 is MTLL, 2 is MTLE, 3 is BG, 4 is TCS, 5 is DU, 6 is vInf

% filtering for the lightning input
inputs.fcut = 300e3;

% magnetic field
inputs.Bmag = 50000e-9;

% cluster and number of nodes
inputs.cluster = 'batchnew';
inputs.numnodes = '8';

% gravity waves! gwave reaches mag (as DN/N0) at maxalt and then stays
% there. lam is horizontal wavelength.

inputs.dogwave = 0;
inputs.gwavemag = 0.5;
inputs.gwavemaxalt = 100e3;
inputs.gwavelam = 20e3;
