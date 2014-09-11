% calculate attenuation due to ionosphere from VLF transmitter runs. Cases
% with no ionosphere provide background; noheating provides attenuation
% without self-consistent heating effects; full has it all.

clear all; close all;

rundir = '/shared/users/ram80/empcodes/runs/transtest/TX_NAA/';

datatype = 'double';

loadconstants;

% parameters
s = get3drunparams(rundir);

pk = max(s.Isource);
% convert source current to power
Txpower = ((2*pi*s.txf0) * s.sourcealt * pk)^2 / (12*pi*e0*vp^3) / 2;

% get probes 
probes = getProbeFields(rundir,s,datatype);

Hpprobe = probes.Hmag(9,:);


%% plot it


h1 = figure(1);
set(h1,'position',[500 200 1000 800]);
set(h1,'Renderer','painters');

ax4 = axes;

tvec = (1:s.tsteps)*s.dt;

plot(ax4,tvec*1e3,Hpprobe*1e6);
xlabel(ax4,'time (ms)');
ylabel(ax4,'H magnitude (uT)');


%% calculate theoretical field along ground from 1 MW transmitter

% assume r = 100 km

dist = 100e3;
Htheo = sqrt(3/4 * e0*vp/pi * Txpower) * 1/dist;

tstart = find(tvec > 0.5e-3,1,'first');
Hpk = max(Hpprobe(tstart:end));

fprintf('Transmitter power = %.1f MW\n',Txpower/1e6);
fprintf('Expected H value at 100 km = %.1f uT\n',Htheo*1e6);
fprintf('Measured H value at 100 km = %.1f uT\n',Hpk*1e6);