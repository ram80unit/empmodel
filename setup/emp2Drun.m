%% setup file: change the parameters as seen fit, then run this script 
% to create the inputs.dat file that will be read by emp2d. 

function [in,jobid] = emp2Drun(in)

loadconstants;

% number of time steps to run
in.maxdist = sqrt(in.range^2 + in.maxalt^2);
in.tsteps = floor(1.1 * in.maxdist / vp / in.dt);

% ground conductivity and epsilon
[in.gsigma,in.gepsilon] = getGroundParams(in);

% option to override conductivity, if in.groundSigma is defined

if isfield(in,'groundSigma'),
  in.gsigma = in.groundSigma * ones(size(in.gsigma));
end

% option to override permittivity, if in.groundEpsilon is defined

  if isfield(in,'groundEpsilon'),
  in.gepsilon = in.groundEpsilon * ones(size(in.gepsilon));
end

% need r vector; this won't be save in the output file, but we need it to
% create an ionosphere.
if (in.groundmethod == 1), 
    in.nground = 0; 
end
[in.r,dr] = generateRvector(in);

in.prober = zeros(size(in.probealt));
in.probet = zeros(size(in.probealt));
in.nprobes = length(in.probealt);
for i = 1:in.nprobes,
    in.prober(i) = find((in.r-in.Re) > in.probealt(i),1,'first') - 1;
    in.probet(i) = floor(in.proberange(i)/in.dr1);
end

% index to start nonlinear calculations
in.nonlinearstart = round(in.nonlinearstartaltitude/in.dr1) + in.nground;

% write out some parameters; approximate guess at run time
hsteps = round(in.range/in.drange);
fprintf('Grid is %d x %d, and will run %d time steps\n',length(in.r),hsteps,in.tsteps);

cellsxtimes = length(in.r)*hsteps*in.tsteps;
timefactor = 3.8e-8;  % empirical
fprintf('Should take about %.1f minutes to run (one node)\n',cellsxtimes*timefactor);


%% input current pulse, defined versus altitude and time. 

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
fwrite(fid,in.numfiles,'int');
fwrite(fid,in.planet,'int');
fwrite(fid,in.decfactor,'int');
fwrite(fid,in.nprobes,'int');
fwrite(fid,in.prober,'int');
fwrite(fid,in.probet,'int');
fwrite(fid,in.dogwave,'int');
fwrite(fid,in.gwavemag,'double');
fwrite(fid,in.gwavemaxalt,'double');
fwrite(fid,in.gwavekh,'double');
fwrite(fid,in.nonlinearstart,'int');
fwrite(fid,in.doDFT,'int');
fwrite(fid,length(in.DFTfreqs),'int');
fwrite(fid,in.DFTfreqs,'double');
fwrite(fid,in.read2Dionosphere,'int');
fclose(fid);


%% write ground parameters to their own file

fid = fopen([in.rundir '/ground.dat'],'w');
fwrite(fid,in.gsigma,'double');
fwrite(fid,in.gepsilon,'double');
fclose(fid);


%% magnetic field specified over 2D domain.

% need to know space parameters
in.rr = in.stepalt/in.dr1 + (in.maxalt - in.stepalt)/in.dr2 + 1 + in.nground;
in.thmax = in.range / in.Re;
in.dth = in.drange / in.Re;
in.hh = round(in.thmax / in.dth) + 1;

% vectors for magnetic field; they will be hh size, then will be made into
% matrices inside code.
in.Br = in.Bmag * ones(1,in.hh); % cos(linspace(0,pi/2,hh));
in.Bt = in.Bmag * zeros(1,in.hh); %sin(linspace(0,pi/2,hh));
in.Bp = in.Bmag * zeros(1,in.hh); %sin(linspace(0,pi/2,hh));

fid = fopen([in.rundir '/B0.dat'],'w');
fwrite(fid,in.Br,'double');
fwrite(fid,in.Bt,'double');
fwrite(fid,in.Bp,'double');
fclose(fid);

nd = MSISatmosphere1((in.r-in.Re)/1000);
ndt = nd.total * 1e6;
in.ndt = ndt;

%% ionosphere and atmosphere densities. Run setupAtmosphere and save ne.dat and nd.dat.

if ~in.read2Dionosphere,
    %ne = IRIDaytime1((in.r-in.Re)/1000);
    %ne = IRIionosphere1((in.r-in.Re)/1000);
    %ne = VictorNeProfile((in.r-in.Re)/1000,1);

    in.ne = YukiIonosphere((in.r-in.Re)/1000,in.beta,in.hprime);
else
    % 2D ionosphere
    in = create2Dionosphere(in);
end


% fix up ne, to continue down to 0 km alt, at which point it will be 1e-5
% electron/m^3

if in.planet == 1,
    in.ne = VenusIonosphere1((in.r-in.Re)/1000,2);
    ndv = VenusAtmosphere((in.r-in.Re)/1000);
    ndt = ndv.total * 1e6;
end

% just for kicks, calculate the maximum dt we could get away with, assuming
% vp dt / ds = sqrt(1 - (wp*dt/2)^2)

ds = 1/(sqrt(1/in.dr2^2 + 1/(in.Re*in.dth)^2));
wpmax = sqrt(max(in.ne(:))*QE^2/ME/e0);
dtmax = 1/sqrt((vp/ds)^2 + (wpmax/2)^2);

fprintf('FYI, maximum usable dt is %.3g us\n',dtmax*1e6);

% ions: same as electrons, except when ne < 100 cm^-3

ni = in.ne;
ni(ni < 100 * 1e6) = 100 * 1e6;

% okay, write them both out to files

fid = fopen([in.rundir '/ne.dat'],'w');
fwrite(fid,in.ne','double');
fclose(fid);

if in.read2Dionosphere && strcmp(in.iono2Dmethod,'eclipse'),
    fid = fopen([in.rundir '/nu.dat'],'w');
    fwrite(fid,in.nu','double');
    fclose(fid);
end

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

% this part is slow; don't do it if rates already exists. Just copy from last run.
if ~isfield(in,'rates'),
    in.rates = getNonlinearRates((in.r-in.Re)/1000,nd,in.nground);
end

rates = in.rates;

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

if in.doplots,
    
h1 = figure(1);
set(h1,'position',[100 100 600 600]);
ax1 = subplot(221);
plot(ax1,log10(in.ne(:,1)),(in.r-in.Re)/1000);
hold(ax1,'on');
%plot(ax1,log10(nec),(r-in.Re)/1000,'.');
plot(ax1,log10(ni),(in.r-in.Re)/1000,'r');
legend(ax1,'electron density','+ion density');

end

% collision frequency for electrons
if in.planet == 0,
    mue = 1.4856 * ndt(in.nground+1) ./ ndt;
else
    mue = 0.0018 * ndt(in.nground+1) ./ ndt;
end
nue = (QE / ME) ./ mue;
nui = nue / 100;

if in.doplots,
    
ax2 = subplot(222);
plot(ax2,log10(nue),(in.r-in.Re)/1000);
hold(ax2,'on');
plot(ax2,log10(nui),(in.r-in.Re)/1000,'r');
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
plot(ax4,nd.temp,(in.r-in.Re)/1000);
xlabel(ax4,'Ambient Temperature');

drawnow;

end

save([in.rundir '/inputs.mat'],'in');

shfile = writeshfile(in);

%% run the simulation

if (in.submitjob),
    
    % run command
    
    system(['cp ' in.exedir in.exefile ' ' in.rundir]);
    
    % for simplicity, cd into run directory, run it, then return to pwd
    
    thisdir = pwd;
    cd(in.rundir);
    
    if strcmp(in.cluster,'local'),
        submitstr = ['sh ' pbsfile ' &'];
    else
        %submitstr = ['qsub -q ' in.cluster ' -d ' in.rundir ' -l nodes=1:ppn=' in.numnodes ' -l walltime=72:00:00 ' ...
        %    pbsfile];
        submitstr = ['sbatch ' shfile];
    end
    
    [~,jobname] = system(submitstr);
    jobid = strtrim(jobname);
    
    fprintf('Job %s submitted!\n',jobid);
    
    cd(thisdir);
    
else
    
    jobid = '';
    fprintf('Job not submitted\n');
    
end
