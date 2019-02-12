clear all; close all;

loadconstants;

% master directory for set of runs
toprundir = '/projects/roma8490/emp/runs/beta0.35/';

% variable for batch of runs. name must match an input!
var1.name = 'hprime';
var1.values = [74 76];

inputs.exefile = 'emp2d';
inputs.exedir = '/projects/roma8490/emp/emp2/';

inputs.submitjob = 1;  % set to zero just to test setup
inputs.savefields = [1 1 1 0 0 0]; % set to zero if you don't need the large output fields
inputs.Re = 6370000;

inputs.dopml_top = 1;
inputs.dopml_wall = 1;
inputs.doionosphere = 1;
inputs.doioniz = 0;
inputs.doelve = 0;
inputs.dodetach = 0;
inputs.numfiles = 50;
inputs.maxalt = 110e3;
inputs.groundmethod = 1; % 0 for PEC, 1 for SIBC, 2 for real

inputs.dr0 = 100;
inputs.nground = 5;
inputs.dr1 = 500;
inputs.dr2 = 250;
inputs.drange = 500;
inputs.stepalt = 70e3;
inputs.nonlinearstartaltitude = 50e3;

inputs.Trlat = 48.204;
inputs.Trlon = -121.917;
inputs.range = 3500e3;
inputs.az = 104.19;

% this is NAA towards Huntsville
%inputs.Trlat = 44.6464;
%inputs.Trlon = -67.2811;
%inputs.range = 2500e3;
%inputs.az = 242.781;

inputs.dt = 1e-7;
inputs.sig = 0;
inputs.sigm = 0;

inputs.camdist = 300e3;
inputs.camalt = 0;
inputs.camelev = elevToElve(inputs.camdist,85e3);
inputs.camfov = [36 18];  %% left-right, up-down
inputs.numpixels = [128 64];   %% camera pixels
inputs.cameratype = 'camera';
inputs.elvesteps = 1000;        % number of elve output time steps
inputs.dotransmitter = 0;
inputs.txf0 = 6e3;
inputs.doDFT = 1;
inputs.DFTfreqs = [5e3 10e3:500:20e3 25e3 30e3];
inputs.planet = 0;  

inputs.decfactor = 4;

inputs.proberange = [1000 2000] * 1e3;
inputs.probealt = [0 0]; %zeros(size(inputs.proberange));

inputs.lightningtype = 0;        % 0 = CG, 1 = my IC, or 2 = Caitano da Silva CID
inputs.I0 = 10e3; % peak current in amperes:
inputs.Ic = 0e3;   % continuing current!
inputs.sourcealt = 5e3; % note that source alt is channel length for CG, altitude of IC
inputs.chlength = 3e3;  % channel length for IC, ignored by CG
inputs.taur = 20e-6;
inputs.tauf = 50e-6;
inputs.rsspeed = -0.75*vp; % for CG, downwards by default; if negative, goes upwards
inputs.decaytype = 1;
% choices: 0 = TL, 1 = MTLL, 2 = MTLE, 3 = BG, 4 = TCS, 5 = DU, 6 = dummy

inputs.fcut = 60e3;

inputs.Bmag = 50000e-9;

inputs.beta = 0.35;
inputs.hprime = 74;

inputs.cluster = 'janus';
inputs.numnodes = '12';
inputs.walltime = '72:00:00';

inputs.dogwave = 0;
inputs.gwavemag = 0.5;
inputs.gwavemaxalt = 100e3;
inputs.gwavekh = 2*pi/20e3;


% submit jobs

for m = 1:length(var1.values),

    % change variables as requested
          evalstr = ['inputs.' var1.name ' = ' num2str(var1.values(m)) ';'];
eval(evalstr);

inputs.runname = [var1.name '_' sprintf('%.2g',var1.values(m))];
inputs.runname(strfind(inputs.runname,'+')) = '';
inputs.rundir = [toprundir inputs.runname];

    % launch job
    [in,jobid] = emp2Drun(inputs);

drawnow;

end
