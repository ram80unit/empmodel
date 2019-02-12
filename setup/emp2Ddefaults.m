% default inputs for emp 2D jobs

inputs.doplots = 0;

inputs.exefile = 'emp2d';
inputs.exedir = '/projects/wexu6668/emp/emp2/';

inputs.submitjob = 1;  % set to zero just to test setup
inputs.savefields = [0 0 0 0 0 0]; % set to zero if you don't need the large output fields
                                   % six files: E, J, H, K, D, O
                                   % K is Eeff, Ek, heat, S, and Te; 
                                   % D is ne, nO-, nu;
                                   % O is optics for seven emissions
% what planet?
inputs.Re = 6370000;

% do you want to use the PML boundary?
inputs.dopml_top = 1;
inputs.dopml_wall = 1;
% do you want to include the ionosphere?
inputs.doionosphere = 1;
% do you want to calculate ionosphere changes (Ne, etc)?
inputs.doioniz = 0;
% do you want to integrate and spit out the elve movie?
inputs.doelve = 0;
% do you want to do associative detachment?
inputs.dodetach = 0;
% number of times to write to output arrays - evenly distributed
inputs.numfiles = 10;
% highest altitude to consider. Notice everything in meters!
inputs.maxalt = 110e3;
% perfectly conducting ground (0), SIBC (1) or real ground (2)? 
inputs.groundmethod = 1;

% dr below surface (ground)
inputs.dr0 = 100;
%number of ground cells
inputs.nground = 5;
% initial dr from ground up
inputs.dr1 = 1000;
% dr at higher altitudes
inputs.dr2 = 500;
% delta in range direction (used to calculate dtheta)
inputs.drange = 1000;
% altitude at which to change to smaller dr. If they are the same, it
% doesn't matter.
inputs.stepalt = 70e3;
% altitude to start nonlinear calculations: avoids heating / ionization
% around source
inputs.nonlinearstartaltitude = 50e3;

% lightning location, azimuth sof interest, and range
    % this is NLK towards Huntsville, AL
%inputs.Trlat = 48.204;
%inputs.Trlon = -121.917;
%inputs.range = 3500e3;
%inputs.az = 104.19;

% this is NAA towards Huntsville
inputs.Trlat = 44.6464;
inputs.Trlon = -67.2811;
inputs.range = 2000e3;
inputs.az = 242.781;

% this is NAU towards Bangor, Maine
%    inputs.Trlat = 18.3987;
%inputs.Trlon = -67.1774;
%inputs.range = 3500e3;   % receiver is at 2957.1 km
%inputs.az = 357.4131;


% time step: use 1e-7 for D-region (good up to 150 km, grids >= 100 m)
inputs.dt = 1e-7;

% conductivities: make them nonzero if you want
inputs.sig = 0;
inputs.sigm = 0;

% camera location, distance from source along ground and altitude
inputs.camdist = 300e3;
inputs.camalt = 0;
inputs.camelev = elevToElve(inputs.camdist,85e3);
inputs.camfov = [36 18];  %% left-right, up-down
inputs.numpixels = [128 64];   %% camera pixels
inputs.cameratype = 'auger';
inputs.elvesteps = 10;        % number of elve output time steps
% transmitter?
inputs.dotransmitter = 0;
inputs.txf0 = 6e3;
inputs.doDFT = 1;
inputs.DFTfreqs = [5e3];
% planet 0 = earth; 1 = venus; 2 = saturn
inputs.planet = 0;  

%decimate outputs before writing. Useful for 100m resolution. Default = 1
inputs.decfactor = 4;
inputs.proberange = [2000] * 1e3;
inputs.probealt = [0]; %zeros(size(inputs.proberange));

% lightning inputs
inputs.lightningtype = 0;        % 0 = CG, 1 = my IC, or 2 = Caitano da Silva CID
inputs.I0 = 10e3; % peak current in amperes:
inputs.Ic = 0e3;   % continuing current!
inputs.sourcealt = 7e3; % note that source alt is channel length for CG, altitude of IC
inputs.chlength = 3e3;  % channel length for IC, ignored by CG
inputs.taur = 20e-6;
inputs.tauf = 50e-6;
inputs.rsspeed = -0.75*vp; % for CG, downwards by default; if negative, goes upwards
inputs.decaytype = 1;
% choices: 0 = TL, 1 = MTLL, 2 = MTLE, 3 = BG, 4 = TCS, 5 = DU, 6 = dummy

% filtering for the lightning input
inputs.fcut = 60e3;

% magnetic field
inputs.Bmag = 50000e-9;

% ionosphere profile
inputs.beta = 0.5;
inputs.hprime = 75;

% 2D ionosphere? specified vs altitude and distance.
inputs.read2Dionosphere = 0;
inputs.iono2Dmethod = 'eclipse';
inputs.iono2Dfile = '/projects/wexu6668/emp/eclipse/CO.mat'
inputs.ionopertdist = 2000e3;

% cluster and number of nodes
inputs.cluster = 'normal';
inputs.partition = 'shas';
inputs.numnodes = '12';
inputs.walltime = '12:00:00';

% gravity waves! gwave reaches mag (as DN/N0) at maxalt and then stays
% there. kh is 2*pi/horizontal wavelength.
inputs.dogwave = 0;
inputs.gwavemag = 0.5;
inputs.gwavemaxalt = 100e3;
inputs.gwavekh = 2*pi/20e3;
