% default inputs for emp 2D jobs

inputs.exefile = 'emp3d';
inputs.exedir = '/shared/users/ram80/empcodes/emp3/';

inputs.submitjob = 1;  % set to zero just to test setup
inputs.savefields = [1 1 1 1 1 1]; % set to zero if you don't need the large output fields

% what planet?
inputs.Re = 6370000;
c = 3e8;
inputs.dovenus = 0;

% do you want to use the PML boundary?
inputs.dopml_top = 1;
inputs.dopml_wall = 1;
% do you want to incl0ude the ionosphere?
inputs.doionosphere = 1;
% do you want to calculate ionosphere changes (Ne, etc)?
inputs.doioniz = 1;
% do you want to integrate and spit out the elve movie?
inputs.doelve = 1;
% do you want to do associative detachment?
inputs.dodetach = 1;
% number of times to write to output arrays - evenly distributed
inputs.numfiles = 30;
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

% lightning location, azimuth of interest, and range
% GW elve:
inputs.Trlat = 44.66;
inputs.Trlon = -102.89;
inputs.range = 100e3;
inputs.az = 298;

% transmitter:

% time step: use 1e-7 for D-region (good up to 150 km, grids >= 100 m)
inputs.dt = 1e-7;

% conductivities: make them nonzero if you want
inputs.sig = 0;
inputs.sigm = 0;

% camera location, distance from source along ground and altitude
inputs.camdist = 474e3;
inputs.camalt = 0;
inputs.camelev = 10;
inputs.camfov = [36 18];  %% left-right, up-down
inputs.numpixels = [128 64];
inputs.cameratype = 'camera';
inputs.elvesteps = 1000;

% transmitter?
inputs.dotransmitter = 0;
inputs.txf0 = 20e3;

%decimate outputs before writing. Useful for 100m resolution. Default = 1.
%Doesn't do anything in 3D code.
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
    inputs.probealt = 100e3 * ones(21,1);
    inputs.proberanget = [0 20 40 60 80 100 -20 -40 -60 -80 -100  0  0  0  0   0   0   0   0   0    0]*1e3;
    inputs.proberangep = [0  0  0  0  0   0   0   0   0   0    0 20 40 60 80 100 -20 -40 -60 -80 -100]*1e3;
end

% lightning inputs
inputs.lightningtype = 0;        % 0 = CG, 1 = IC, or 2 = CID
inputs.sourcedirection = 0;      % 0 is vertical, 1 is horizontal, only applies if IC or CID
inputs.I0 = 100e3; % peak current in amperes:
inputs.Ic = 0; %2e3;   % continuing current!
inputs.sourcealt = 8e3;
inputs.chlength = 2e3;  % channel length for IC, ignored by CG
inputs.taur = 10e-6;
inputs.tauf = 50e-6;
inputs.rsspeed = 0.7*c;
inputs.decaytype = 1;
% 0 is TL, 1 is MTLL, 2 is MTLE, 3 is BG, 4 is TCS, 5 is DU, 6 is vInf

% filtering for the lightning input
inputs.fcut = 300e3;

% magnetic field: r, theta, phi components
%inputs.Bmag = [50000e-9 0 0];
%inputs.Bmag = [-49740e-9 17145e-9 -37759e-9];
inputs.Bdip = NaN;      % temporary, otherwise Bmag gets overwritten later
[Bx, By, Bz] = igrf(today, inputs.Trlat, inputs.Trlon, 90, 'geodetic');
%inputs.Bmag = [-Bz Bx By] * 1e-9;
% 45-degree field with 50,000 nT amplitude
inputs.Bmag = [-35355 0 35355] * 1e-9;

% cluster and number of nodes
inputs.cluster = 'batchnew';
inputs.numnodes = '12';

% gravity waves! gwave reaches mag (as DN/N0) at maxalt and then stays
% there. kh and kp are horizontal k (2-pi/wavelength).

inputs.dogwave = 0;
inputs.gwavemag = 0.2;
inputs.gwavemaxalt = 90e3;
inputs.gwavekr = 2*pi/20e3;
inputs.gwavekh = 2*pi/60e3;
inputs.gwavekp = 2*pi/500e3;

% modified IRI profile
inputs.doIRI = 1;
% relative factor to multiply my default IRI profile
inputs.IRI = 1;
