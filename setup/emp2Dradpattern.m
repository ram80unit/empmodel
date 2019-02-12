%% setup file: change the parameters as seen fit, then run this script 
% to create the inputs.dat file that will be read by emp2d. 

% this version: determine radiation pattern, no ionosphere. test different
% return stroke models.

clear all; close all;

loadconstants;

runname = 'tcstest';
rundir = ['/shared/users/ram80/empcodes/runs/radpattern/' runname '/'];
exefile = 'emp2d';
exedir = '/shared/users/ram80/empcodes/emp2/';

submitjob = 1;  % set to zero just to test setup


% do you want to use the PML boundary?
dopml_top = 1;
dopml_wall = 1;
% do you want to include the ionosphere?
doionosphere = 0;
% do you want to calculate ionosphere changes (Ne, etc)?
doioniz = 0;
% do you want to integrate and spit out the elve movie?
doelve = 0;
% do you want to do associative detachment?
dodetach = 0;
% number of times to write to output arrays - evenly distributed
numfiles = 60;
% highest altitude to consider. Notice everything in meters!
maxalt = 100e3;
% perfectly conducting ground, or real ground? 
pecground = 1;

% dr below surface (ground)
dr0 = 10;
%number of ground cells
nground = 5;
% initial dr from ground up
dr1 = 200;
% dr at higher altitudes
dr2 = 200;
% altitude at which to change to smaller dr. If they are the same, it
% doesn't matter.
stepalt = 70e3;

% lightning location, azimuth of interest, and range
Trlat = 40;
Trlon = -100;
range = 100e3;
az = 298;

% time step: use 1e-7 for D-region (good up to 150 km, grids >= 100 m)
dt = 1e-7;

% number of time steps to run
maxdist = sqrt(range^2 + maxalt^2);
tsteps = floor(1.5 * maxdist / vp / dt);

% conductivities: make them nonzero if you want
sig = 0;
sigm = 0;

% ground conductivity and epsilon
[gsigma,gepsilon] = getGroundParams(Trlat,Trlon,az,range,dr1);

% camera location, distance from source along ground and altitude
camdist = 400e3;
camalt = 0;

% transmitter? venus?
dotransmitter = 0;
dovenus = 0;
txf0 = 20e3;

%decimate outputs before writing. Useful for 100m resolution. Default = 1
decfactor = 2;

% set up probe points. Define locations in km, then determine grid values
probedist = 70; 
probeangle = linspace(0,90,31);

probealt = probedist * sind(probeangle);
proberange = probedist * cosd(probeangle);
 
% need r vector; this won't be save in the output file, but we need it to
% create an ionosphere.
[r,dr] = generateRvector(dr0,dr1,dr2,nground,stepalt,maxalt);

prober = zeros(size(probealt));
probet = zeros(size(probealt));
nprobes = length(probealt);
for i = 1:nprobes,
    prober(i) = find((r-RE) > probealt(i)*1e3,1,'first') - 1;
    probet(i) = floor(proberange(i)*1e3/dr1);
end

% write out some parameters; approximate guess at run time
fprintf('Grid is %d x %d, and will run %d time steps\n',length(r),round(range/dr1),tsteps);

cellsxtimes = length(r)*round(range/dr1)*tsteps;
timefactor = 3.8e-8;  % empirical
fprintf('Should take about %.1f minutes to run (one node)\n',cellsxtimes*timefactor);


%% input current pulse (vs time) for lightning strike. This can be anything, as long as it has tsteps values. Here, we implement some filtering to keep dr < 10*lambda_max. 

Iin = zeros(tsteps,1);

% Source waveform.
% for NPM, use I0 = 160 A, L = 2 km, alt = 1 km

I0 = 100e3; % peak current in amperes:
Ic = 0; %2e3;   % continuing current!
sourcealt = 20e3;
taur = 10e-6;
tauf = 30e-6;
rsspeed = 2e8;
% current source type: 0 = TL, no decay; 1 = MTLL, 2 = MTLE
% 3 = BG, 4 = TCS, 5 = DU
% 6 is constant-current dipole (v = infty)
sourcetype = 4;


% okay, set up current waveform: linear rise, exp decay
for t = 1:tsteps,
    if (t * dt < taur),
        Iin(t) = I0 * t * dt / taur;
    else
        Iin(t) = I0 * exp(-(t*dt-taur)^2/tauf^2) + Ic*(1 - exp(-(t*dt-taur)^2/tauf^2));
    end
end

% okay, filter it

fcut = 300e3;
[b,a] = butter(2,fcut*dt/2);

Iin2 = filter(b,a,Iin);
Iin2 = Iin2 / max(Iin2) * I0;   % rescaled to peak current

% fix delays and stuff
Iin2 = [zeros(2*round(sourcealt/rsspeed/dt),1); Iin2];
Iin2(Iin2 < 1e-6) = 1e-6;

tsteps = length(Iin2);

% write everything to a file, to be read by program

if ~exist(rundir,'dir'),
    mkdir(rundir);
end

fid = fopen([rundir '/inputs.dat'],'w');
fwrite(fid,dopml_top,'int');
fwrite(fid,dopml_wall,'int');
fwrite(fid,doionosphere,'int');
fwrite(fid,doioniz,'int');
fwrite(fid,doelve,'int');
fwrite(fid,dodetach,'int');
fwrite(fid,dotransmitter,'int');
fwrite(fid,pecground,'int');
fwrite(fid,maxalt,'double');
fwrite(fid,stepalt,'double');
fwrite(fid,dr0,'double');
fwrite(fid,dr1,'double');
fwrite(fid,dr2,'double');
fwrite(fid,nground,'int');
fwrite(fid,range,'double');
fwrite(fid,dt,'double');
fwrite(fid,tsteps,'int');
fwrite(fid,sig,'double');
fwrite(fid,sigm,'double');
fwrite(fid,camdist,'double');
fwrite(fid,camalt,'double');
fwrite(fid,Iin2,'double');
fwrite(fid,txf0,'double');
fwrite(fid,sourcealt,'double');
fwrite(fid,rsspeed,'double');
fwrite(fid,sourcetype,'int');
fwrite(fid,numfiles,'int');
fwrite(fid,dovenus,'int');
fwrite(fid,decfactor,'int');
fwrite(fid,nprobes,'int');
fwrite(fid,prober,'int');
fwrite(fid,probet,'int');
fclose(fid);


%% write ground parameters to their own file

fid = fopen([rundir 'ground.dat'],'w');
fwrite(fid,gsigma,'double');
fwrite(fid,gepsilon,'double');
fclose(fid);


%% magnetic field specified over 2D domain.

% need to know space parameters
rr = stepalt/dr1 + (maxalt - stepalt)/dr2 + 1 + nground;
thmax = range / RE;
dth = dr1 / RE;
hh = round(thmax / dth) + 1;

% temporarily testing basic homogeneous field
Bmag = 50000e-9;
if dovenus,
    Bmag = 50e-9;
end

% vectors for magnetic field; they will be hh size, then will be made into
% matrices inside code.
Br = Bmag * ones(1,hh); % cos(linspace(0,pi/2,hh));
Bt = Bmag * zeros(1,hh); %sin(linspace(0,pi/2,hh));
Bp = Bmag * zeros(1,hh); %sin(linspace(0,pi/2,hh));

fid = fopen([rundir 'B0.dat'],'w');
fwrite(fid,Br,'double');
fwrite(fid,Bt,'double');
fwrite(fid,Bp,'double');
fclose(fid);


%% ionosphere and atmosphere densities. Run setupAtmosphere and save ne.dat and nd.dat.

beta = 0.8;
hk = 82;

%ne = IRIDaytime1((r-RE)/1000);
ne = IRIionosphere1((r-RE)/1000);
%ne = VictorNeProfile((r-RE)/1000,1);
%ne = YukiIonosphere((r-RE)/1000,beta,hk);
%nec = YukiIonosphere((r-RE)/1000,beta,85);
nd = MSISatmosphere1((r-RE)/1000);
ndt = nd.total * 1e6;

if dovenus,
    ne = VenusIonosphere1((r-RE)/1000,2);
    ndv = VenusAtmosphere((r-RE)/1000);
    ndt = ndv.total * 1e6;
end

% just for kicks, calculate the maximum dt we could get away with, assuming
% vp dt / ds = sqrt(1 - (wp*dt/2)^2)

dth = dr1/RE;
ds = 1/(sqrt(1/dr2^2 + 1/(RE*dth)^2));
wpmax = sqrt(max(ne)*QE^2/ME/e0);
dtmax = 1/sqrt((vp/ds)^2 + (wpmax/4)^2);

fprintf('FYI, maximum usable dt is %.3g us\n',dtmax*1e6);

% ions: same as electrons, except when ne < 100 cm^-3

ni = ne;
ni(ni < 100 * 1e6) = 100 * 1e6;

% okay, write them both out to files

fid = fopen([rundir 'ne.dat'],'w');
fwrite(fid,ne,'double');
fclose(fid);

fid = fopen([rundir 'nd.dat'],'w');
fwrite(fid,ndt,'double');
fclose(fid);

fid = fopen([rundir 'ni.dat'],'w');
fwrite(fid,ni,'double');
fclose(fid);

fid = fopen([rundir 'etemp.dat'],'w');
fwrite(fid,nd.temp,'double');
fclose(fid);

% set up rates for ionization, attachment, mobility, and optics. These will
% then be read and interpolated in the code.

rates = getNonlinearRates((r-RE)/1000,nd);

% save to file.

fid = fopen([rundir 'rates.dat'],'w');
fwrite(fid,length(rates.efield),'int');
fwrite(fid,rates.efield,'double');
fwrite(fid,rates.ioniz,'double');
fwrite(fid,rates.attach,'double');
fwrite(fid,rates.mobility,'double');
fwrite(fid,rates.Ored,'double');
fwrite(fid,rates.Ogrn,'double');
fwrite(fid,rates.N21p,'double');
fwrite(fid,rates.N22p,'double');
fwrite(fid,rates.N2p1N,'double');
fwrite(fid,rates.N2pM,'double');
fwrite(fid,rates.O2p1N,'double');
fclose(fid);


%% quick plots

h1 = figure(1);
set(h1,'position',[100 100 600 600]);
ax1 = subplot(221);
plot(ax1,log10(ne),(r-RE)/1000);
hold(ax1,'on');
%plot(ax1,log10(nec),(r-RE)/1000,'.');
plot(ax1,log10(ni),(r-RE)/1000,'r');
legend(ax1,'electron density','+ion density');

% collision frequency for electrons
if ~dovenus,
    mue = 1.4856 * ndt(nground+1) ./ ndt;
else
    mue = 0.0018 * ndt(nground+1) ./ ndt;
end
nue = (QE / ME) ./ mue;
nui = nue / 100;

ax2 = subplot(222);
plot(ax2,log10(nue),(r-RE)/1000);
hold(ax2,'on');
plot(ax2,log10(nui),(r-RE)/1000,'r');
legend(ax2,'electron coll. freq.','ion coll. freq.');

% source

ax3 = subplot(223);
tvec = 0:dt:(length(Iin2)-1)*dt;
plot(ax3,tvec*1e6,Iin2/1000);
set(ax3,'xlim',[0 tsteps*2*dt*1e6]);
xlabel(ax3,'time (us)');
ylabel(ax3,'Current (kA)');

ax4 = subplot(224);
plot(ax4,nd.temp,(r-RE)/1000);
xlabel(ax4,'Ambient Temperature');


%% run the simulation

if (submitjob),
    
    % create pbs file to run simulation
    
    pbsfile = writepbsfile(rundir,runname,exefile);
    
    % run command
    
    system(['cp ' exedir exefile ' ' rundir]);
    
    % for simplicity, cd into run directory, run it, then return to pwd
    
    thisdir = pwd;
    cd(rundir);
    
    submitstr = ['qsub -q batchnew -o ' rundir 'log.txt -d ' rundir ' -l nodes=1:ppn=8 -l walltime=72:00:00 ' ...
        pbsfile];
    
    [tmp,jobname] = system(submitstr);
    jobname = strtrim(jobname);
    
    fprintf('Job %s submitted!\n',jobname);
    
    cd(thisdir);
    
end
