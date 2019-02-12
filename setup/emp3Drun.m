%% setup file: change the parameters as seen fit, then run this script 
% to create the inputs.dat file that will be read by emp2d. 

function [in,jobid] = emp3Drun(in)

loadconstants;

% number of time steps to run
in.maxdist = sqrt(2*(in.range^2) + in.maxalt^2);
in.tsteps = floor(1.1 * in.maxdist / vp / in.dt);

% ground conductivity and epsilon
% in 3D it is just scalar values
[gsigma,gepsilon] = getGroundParams(in.Trlat,in.Trlon,in.az,in.range,in.dr1);
in.gsigma = gsigma(1);
in.gepsilon = gepsilon(1); 

% need r vector; this won't be save in the output file, but we need it to
% create an ionosphere.
if (in.groundmethod == 1), 
    in.nground = 0; 
end
[in.r,dr] = generateRvector(in.dr0,in.dr1,in.dr2,in.nground,in.stepalt,in.maxalt);
hh = 2*floor(in.range/in.drange);
pp = 2*floor(in.range/in.drange);

% set up probe points. Define locations in km, then determine grid values

in.prober = zeros(size(in.probealt));
in.probet = zeros(size(in.probealt));
in.probep = zeros(size(in.probealt));
in.nprobes = length(in.probealt);
for i = 1:in.nprobes,
    in.prober(i) = find((in.r-RE) > in.probealt(i),1,'first') - 1;
    in.probet(i) = floor(hh/2 + in.proberanget(i)/in.drange);
    in.probep(i) = floor(pp/2 + in.proberangep(i)/in.drange);
end


% write out some parameters; approximate guess at run time
hsteps = round(2*in.range/in.drange + 1);
fprintf('Grid is %d x %d x %d, and will run %d time steps\n',length(in.r),hsteps,hsteps,in.tsteps);

% need to recalculate this for 3D
cellsxtimes = length(in.r)*hsteps*hsteps*in.tsteps/str2num(in.numnodes);
timefactor = 6.63e-8;  % empirical
fprintf('Should take about %.1f minutes to run (%s nodes)\n',cellsxtimes*timefactor,in.numnodes);


%% input current pulse, defined vs. altitude and time. 

% set up source file: will be an array of doubles, size (source alt) x tsteps.
in = createEmpSource(in);
nt_source = size(in.source,2);
nalt_source = size(in.source,1);

% code needs to know channel length for ICs, in number of cells. I will
% give it HALF channel length!
channelcells = round(in.chlength/in.dr1/2);


%% write everything to a file, to be read by program

if ~exist(in.rundir,'dir'),
    mkdir(in.rundir);
end

fid = fopen([in.rundir '/source.dat'],'w');
fwrite(fid,nalt_source,'int');
fwrite(fid,nt_source,'int');
fwrite(fid,channelcells,'int');
fwrite(fid,in.source,'double');
fclose(fid);

% regular input file

fid = fopen([in.rundir '/inputs.dat'],'w');
fwrite(fid,in.Re,'double');
fwrite(fid,in.dopml_top,'int');
fwrite(fid,in.dopml_wall,'int');
fwrite(fid,in.doionosphere,'int');
fwrite(fid,in.doioniz,'int');
fwrite(fid,in.doelve,'int');
fwrite(fid,in.dodetach,'int');
fwrite(fid,in.dotransmitter,'int');
fwrite(fid,in.savefields,'int');
fwrite(fid,in.groundmethod,'int');
fwrite(fid,in.maxalt,'double');
fwrite(fid,in.stepalt,'double');
fwrite(fid,in.dr0,'double');
fwrite(fid,in.dr1,'double');
fwrite(fid,in.dr2,'double');
fwrite(fid,in.nground,'int');
fwrite(fid,in.range,'double');
fwrite(fid,in.drange,'double');
fwrite(fid,in.dt,'double');
fwrite(fid,in.tsteps,'int');
fwrite(fid,in.sig,'double');
fwrite(fid,in.sigm,'double');
fwrite(fid,in.camdist,'double');
fwrite(fid,in.camalt,'double');
fwrite(fid,in.elvesteps,'int');
fwrite(fid,in.sourcedirection,'int');
fwrite(fid,in.numfiles,'int');
fwrite(fid,in.dovenus,'int');
fwrite(fid,in.decfactor,'int');
fwrite(fid,in.nprobes,'int');
fwrite(fid,in.prober,'int');
fwrite(fid,in.probet,'int');
fwrite(fid,in.probep,'int');
fwrite(fid,in.dogwave,'int');
fwrite(fid,in.gwavemag,'double');
fwrite(fid,in.gwavemaxalt,'double');
fwrite(fid,in.gwavekr,'double');
fwrite(fid,in.gwavekh,'double');
fwrite(fid,in.gwavekp,'double');
fclose(fid);

if in.dogwave,
    gwavelam = 2*pi/sqrt(in.gwavekh^2 + in.gwavekp^2);
    fprintf('Gravity Wave! Delta amplitude at %.0f km altitude is %.2f, horiz wavelength is %.1f km\n',...
        in.gwavemaxalt/1e3,in.gwavemag,gwavelam/1e3);
end

%% write ground parameters to their own file

fid = fopen([in.rundir '/ground.dat'],'w');
fwrite(fid,in.gsigma,'double');
fwrite(fid,in.gepsilon,'double');
fclose(fid);


%% if do transmitter, print out (and save) the total radiated power, assuming a small dipole

if in.dotransmitter,
    omega = 2*pi*in.txf0;
    in.Prad = (omega * in.sourcealt * max(in.Iin))^2 / (12 * pi * e0 * vp^3) / 2;
    
    fprintf('Transmitter: total radiated power is %.2f MW\n',in.Prad/1e6);
end


%% magnetic field specified over 2D domain.

% need to know space parameters
in.rr = in.stepalt/in.dr1 + (in.maxalt - in.stepalt)/in.dr2 + 1 + in.nground;
in.thmax = in.range / RE;
in.phmax = in.range / RE;
in.dth = in.drange / RE;
in.dph = in.drange / RE;
in.hh = 2 * round(in.thmax / in.dth) + 1;
in.pp = 2 * round(in.phmax / in.dph) + 1;

% changing dip angle requires some modifications from defaults
if ~isnan(in.Bdip),
    Bamp = norm(in.Bmag);
    in.Bmag = [-Bamp*cosd(in.Bdip) 0 Bamp*sind(in.Bdip)];
end

% homogeneous magnetic field in 3D
in.Br = in.Bmag(1);
in.Bt = in.Bmag(2);
in.Bp = in.Bmag(3);

fid = fopen([in.rundir '/B0.dat'],'w');
fwrite(fid,in.Br,'double');
fwrite(fid,in.Bt,'double');
fwrite(fid,in.Bp,'double');
fclose(fid);


%% ionosphere and atmosphere densities. Run setupAtmosphere and save ne.dat and nd.dat.

%beta = 0.9;
%hk = 87;

if in.doIRI,
    ne = IRIionosphere1((in.r-RE)/1000);
    ne = ne * in.IRI;
else
    ne = YukiIonosphere((in.r-RE)/1000,in.beta,in.hk);
end

% some other possibilities:
%ne = IRIDaytime1((r-RE)/1000);
%ne = VictorNeProfile((r-RE)/1000,1);

nd = MSISatmosphere1((in.r-RE)/1000);
ndt = nd.total * 1e6;

% fix up ne, to continue down to 0 km alt, at which point it will be 1e-5
% electron/m^3

if in.dovenus,
    ne = VenusIonosphere1((in.r-RE)/1000,2);
    ndv = VenusAtmosphere((in.r-RE)/1000);
    ndt = ndv.total * 1e6;
end

% just for kicks, calculate the maximum dt we could get away with, assuming
% vp dt / ds = sqrt(1 - (wp*dt/2)^2)

ds = 1/(sqrt(1/in.dr2^2 + 1/(RE*in.dth)^2 + 1/(RE*in.dph)^2));
wpmax = sqrt(max(ne)*QE^2/ME/e0);
dtmax = 1/sqrt((vp/ds)^2 + (wpmax/4)^2);

fprintf('FYI, maximum usable dt is %.3g us\n',dtmax*1e6);

% ions: same as electrons, except when ne < 100 cm^-3

ni = ne;
ni(ni < 100 * 1e6) = 100 * 1e6;

% okay, write them both out to files

fid = fopen([in.rundir '/ne.dat'],'w');
fwrite(fid,ne,'double');
fclose(fid);

fid = fopen([in.rundir '/nd.dat'],'w');
fwrite(fid,ndt,'double');
fclose(fid);

fid = fopen([in.rundir '/ni.dat'],'w');
fwrite(fid,ni,'double');
fclose(fid);

fid = fopen([in.rundir '/etemp.dat'],'w');
fwrite(fid,nd.temp,'double');
fclose(fid);

%% set up rates for ionization, attachment, mobility, and optics. These will
% then be read and interpolated in the code.

rates = getNonlinearRates((in.r-RE)/1000,nd,in.nground);

% save to file.

fid = fopen([in.rundir '/rates.dat'],'w');
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

in.rates = rates;


%% camera parameters

fov.leftright = in.camfov(1);
fov.updown = in.camfov(2);
fov.elevation = in.camelev;
fov.numpixels = in.numpixels;

[az,el,totalpixels] = getCameraPixels(in.cameratype,fov);

fid = fopen([in.rundir '/camera.dat'],'w');
if strcmp(in.cameratype,'camera'),
    fwrite(fid,1,'int');
else
    fwrite(fid,0,'int');
end
fwrite(fid,totalpixels,'int');
fwrite(fid,az,'double');
fwrite(fid,el,'double');
fclose(fid);


%% quick plots

h1 = figure(1);
set(h1,'position',[100 100 600 600]);
ax1 = subplot(221);
plot(ax1,log10(ne),(in.r-RE)/1000);
hold(ax1,'on');
%plot(ax1,log10(nec),(r-RE)/1000,'.');
plot(ax1,log10(ni),(in.r-RE)/1000,'r');
legend(ax1,'electron density','+ion density');

% collision frequency for electrons
if ~in.dovenus,
    mue = 1.4856 * ndt(in.nground+1) ./ ndt;
else
    mue = 0.0018 * ndt(in.nground+1) ./ ndt;
end
nue = (QE / ME) ./ mue;
nui = nue / 100;

ax2 = subplot(222);
plot(ax2,log10(nue),(in.r-RE)/1000);
hold(ax2,'on');
plot(ax2,log10(nui),(in.r-RE)/1000,'r');
legend(ax2,'electron coll. freq.','ion coll. freq.');

% source

ax3 = subplot(223);
tvec = in.dt*(0:1:(nt_source-1));
hvec = in.dr1*(0:1:(nalt_source-1));
imagesc(tvec*1e6,hvec/1e3,in.source,'parent',ax3);
axis(ax3,'xy');
xlabel(ax3,'time (us)');
ylabel(ax3,'Altitude (km)');

ax4 = subplot(224);
plot(ax4,nd.temp,(in.r-RE)/1000);
xlabel(ax4,'Ambient Temperature');

drawnow;

save([in.rundir '/inputs.mat'],'in');

%% run the simulation

if (in.submitjob),
    
    % create pbs file to run simulation
    
    pbsfile = writepbsfile(in.rundir,in.runname,in.exefile);
    
    % run command
    
    system(['cp ' in.exedir in.exefile ' ' in.rundir]);
    
    % for simplicity, cd into run directory, run it, then return to pwd
    
    thisdir = pwd;
    cd(in.rundir);
    
    submitstr = ['qsub -q ' in.cluster ' -d ' in.rundir ' -l nodes=1:ppn=' in.numnodes ' -l walltime=72:00:00 ' ...
        pbsfile];
    
    [~,jobname] = system(submitstr);
    jobid = strtrim(jobname);
    
    fprintf('Job %s submitted!\n',jobid);
    
    cd(thisdir);
    
else
    
    jobid = '';
    fprintf('Job not submitted\n');
    
end
