%% setup file: change the parameters as seen fit, then run this script 
% to create the inputs.dat file that will be read by emp2d. 

clear all; close all;

loadconstants;
topdir = '/shared/users/ram80/empcodes/emp2/latest/';

% 2D or 3D?
dims = 2;

% do you want to use the PML boundary?
dopml_top = 1;
dopml_wall = 1;
% do you want to calculate ionosphere changes (Ne, etc)?
doioniz = 1;
% do you want to integrate and spit out the elve movie?
doelve = 1;
% do you want to do associative detachment?
dodetach = 1;
% number of files to write out
numfiles = 60;
% highest altitude to consider. Notice everything in meters!
maxalt = 110e3;

% initial dr from ground up
dr1 = 1000;
% dr at higher altitudes
dr2 = 200;
% altitude at which to change to smaller dr. If they are the same, it
% doesn't matter.
stepalt = 70e3;

% lightning location, azimuth of interest, and range
Trlat = 35;
Trlon = -100;
range = 120e3;
az = 0;

% time step: use 1e-7 for D-region (good up to 150 km)
dt = 1e-7;

% number of time steps to run
maxdist = sqrt(range^2 + maxalt^2);
tsteps = floor(1.3 * maxdist / vp / dt);

% conductivities: make them nonzero if you want
sig = 0;
sigm = 0;

% ground conductivity and epsilon
[gsigma,gepsilon] = getGroundParams(Trlat,Trlon,az,range,dr1);

% camera location, distance from source along ground and altitude
camdist = 400e3;
camalt = 0;


%% input current pulse (vs time) for lightning strike. This can be anything, as long as it has tsteps values. Here, we implement some filtering to keep dr < 10*lambda_max. 

Iin = zeros(tsteps,1);

% Source waveform.
% for NPM, use I0 = 160 A, L = 2 km, alt = 1 km

I0 = 150e3; % peak current in amperes:
Ic = 0; %2e3;   % continuing current!
sourcealt = 10e3;
taur = 20e-6;
tauf = 60e-6;
rsspeed = 2e8;

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

if ~exist(topdir,'dir'),
    mkdir(topdir);
end

fid = fopen([topdir '/inputs.dat'],'w');
fwrite(fid,dopml_top,'int');
fwrite(fid,dopml_wall,'int');
fwrite(fid,doioniz,'int');
fwrite(fid,doelve,'int');
fwrite(fid,dodetach,'int');
fwrite(fid,maxalt,'double');
fwrite(fid,stepalt,'double');
fwrite(fid,dr1,'double');
fwrite(fid,dr2,'double');
fwrite(fid,range,'double');
fwrite(fid,dt,'double');
fwrite(fid,tsteps,'int');
fwrite(fid,sig,'double');
fwrite(fid,sigm,'double');
fwrite(fid,Bmag,'double');
fwrite(fid,Bvec,'double');
fwrite(fid,camdist,'double');
fwrite(fid,camalt,'double');
fwrite(fid,Iin2,'double');
fwrite(fid,sourcealt,'double');
fwrite(fid,rsspeed,'double');
fwrite(fid,numfiles,'int');
fclose(fid);


%% write ground parameters to their own file

fid = fopen([topdir 'ground.dat'],'w');
% 2d will use complete ground; 3d will just use value at source
if dims < 3,
    fwrite(fid,gsigma,'double');
    fwrite(fid,gepsilon,'double');
else
    fwrite(fid,gsigma(1),'double');
    fwrite(fid,gepsilon(1),'double');
end
fclose(fid);


%% ionosphere and atmosphere densities. Run setupAtmosphere and save ne.dat and nd.dat.

% need r vector; this won't be save in teh output file, but we need it to
% create an ionosphere.

[r,dr,dr0,nground] = generateRvector(dr1,dr2,stepalt,maxalt);

beta = 0.8;

%ne = IRIDaytime1((r-RE)/1000);
ne = IRIionosphere1((r-RE)/1000);
%ne = VictorNeProfile((r-RE)/1000,1);
%ne = YukiIonosphere((r-RE)/1000,beta,hk);
%nec = YukiIonosphere((r-RE)/1000,beta,85);
nd = MSISatmosphere1((r-RE)/1000);
ndt = nd.total * 1e6;

% fix up ne, to continue down to 0 km alt, at which point it will be 1e-5
% electron/m^3

fi = find(ne > 0,1,'first');
logne = interp1([1 fi],[-5 log10(ne(fi))],1:fi,'linear');
ne(1:fi) = 10.^(logne);

% just for kicks, calculate the maximum dt we could get away with, assuming
% vp dt / ds = sqrt(1 - (wp*dt/2)^2)

dth = dr1/RE;

ds = 1/(sqrt(1/dr2^2 + 1/(RE*dth)^2));
wpmax = sqrt(max(ne)*QE^2/ME/e0);

dtmax = 1/sqrt((vp/ds)^2 + (wpmax/4)^2);

fprintf('FYI, maximum usable dt is %.3g us\n',dtmax*1e6);

% ions: same as electrons, except when ne < 100 cm^-3

ni = ne;
ni(ni < 10 * 1e6) = 10 * 1e6;

% okay, write them both out to files

fid = fopen([topdir 'ne.dat'],'w');
fwrite(fid,ne,'double');
fclose(fid);

fid = fopen([topdir 'nd.dat'],'w');
fwrite(fid,ndt,'double');
fclose(fid);

fid = fopen([topdir 'ni.dat'],'w');
fwrite(fid,ni,'double');
fclose(fid);

fid = fopen([topdir 'etemp.dat'],'w');
fwrite(fid,nd.temp,'double');
fclose(fid);

% set up rates for ionization, attachment, mobility, and optics. These will
% then be read and interpolated in the code.

rates = getNonlinearRates((r-RE)/1000,nd);

% save to file.
% write order: int number of e-field steps; float e-fields;
% float arrays

fid = fopen([topdir 'rates.dat'],'w');
fwrite(fid,length(rates.efield),'int');
fwrite(fid,rates.efield,'double');
fwrite(fid,rates.ioniz,'double');
fwrite(fid,rates.attach,'double');
fwrite(fid,rates.mobility,'double');
fwrite(fid,rates.Ored,'double');
fwrite(fid,rates.Ogrn,'double');
fwrite(fid,rates.N21p,'double');
fwrite(fid,rates.N22p,'double');
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

mue = 1.4856 * ndt(nground+1) ./ ndt;
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
