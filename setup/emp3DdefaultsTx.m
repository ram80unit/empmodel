% default inputs for emp 2D jobs

inputs.exefile = 'emp3d';
inputs.exedir = '/shared/users/ram80/empcodes/emp3/';

inputs.submitjob = 1;  % set to zero just to test setup
inputs.savefields = 1; % set to zero if you don't need the large output fields

% what planet?
inputs.Re = 6370000;
c = 3e8;

% do you want to use the PML boundary?
inputs.dopml_top = 1;
inputs.dopml_wall = 1;
% do you want to incl0ude the ionosphere?
inputs.doionosphere = 1;
% do you want to calculate ionosphere changes (Ne, etc)?
inputs.doioniz = 1;
% do you want to integrate and spit out the elve movie?
inputs.doelve = 0;
% do you want to do associative detachment?
inputs.dodetach = 1;
% number of times to write to output arrays - evenly distributed
inputs.numfiles = 30;
% highest altitude to consider. Notice everything in meters!
inputs.maxalt = 110e3;
% perfectly conducting ground, or real ground? 
inputs.pecground = 0;

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

% location, azimuth of interest, and range
% transmitter: NPM
inputs.Trlat = 21.420383;
inputs.Trlon = -158.153986;
inputs.range = 200e3;
inputs.az = 180;

% time step: use 1e-7 for D-region (good up to 150 km, grids >= 100 m)
inputs.dt = 2e-7;

% conductivities: make them nonzero if you want
inputs.sig = 0;
inputs.sigm = 0;

% camera location, distance from source along ground and altitude
inputs.camdist = 400e3;
inputs.camalt = 0;
inputs.camelev = 10;
inputs.camfov = [36 18];  %% left-right, up-down
inputs.numpixels = [128 64];
inputs.cameratype = 'camera';
inputs.elvesteps = 1000;
% transmitter? venus?
inputs.dotransmitter = 1;
inputs.dovenus = 0;
inputs.txf0 = 21.4e3;

%decimate outputs before writing. Useful for 100m resolution. Default = 1
inputs.decfactor = 1;

% probe distance for radiation pattern
% for radiation pattern do the following;
if (0),
    inputs.probedist = 90e3;
    probeangle = linspace(0,90,10);
    inputs.probealt = inputs.probedist * sind(probeangle);
    inputs.proberanget = inputs.probedist * cosd(probeangle);
    inputs.proberangep = zeros(size(probeangle));
else
    % set up hard-coded probe alts and ranges
    inputs.probealt = [100e3 * ones(8,1); 0];
    inputs.proberanget = [50 0 -50 0 100 0 -100 0 100]*1e3;
    inputs.proberangep = [0 -50 0 50 0 -100 0 100 0]*1e3;
end

% lightning inputs
inputs.I0 = 200; % peak current in amperes:
inputs.Ic = 0; %2e3;   % continuing current!
inputs.sourcealt = 4e3;
inputs.taur = 10e-6;
inputs.tauf = 50e-6;
inputs.rsspeed = 0.75*c;
inputs.decaytype = 1;
% 0 is TL, 1 is MTLL, 2 is MTLE, 3 is BG, 4 is TCS, 5 is DU, 6 is vInf

% filtering for the lightning input
inputs.fcut = 300e3;

% magnetic field: r, theta, phi components
% get B-field components from igrf
% Bx is northward, By is eastward, Bz is downward
inputs.Bdip = NaN;      % temporary, otherwise Bmag gets overwritten later
[Bx, By, Bz] = igrf(today, inputs.Trlat, inputs.Trlon, 90, 'geodetic');
inputs.Bmag = [-Bz Bx By] * 1e-9;

% cluster and number of nodes
inputs.cluster = 'batch';
inputs.numnodes = '8';

% gravity waves! gwave reaches mag (as DN/N0) at maxalt and then stays
% there. kh and kp are horizontal k (2-pi/wavelength).

inputs.dogwave = 0;
inputs.gwavemag = 0.5;
inputs.gwavemaxalt = inputs.maxalt;
inputs.gwavekr = 2*pi/20e3;
inputs.gwavekh = 2*pi/60e3;
inputs.gwavekp = 2*pi/120e3;

% modified IRI profile
inputs.doIRI = 1;
% relative factor to multiply my default IRI profile
inputs.IRI = 1;
