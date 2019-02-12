%% setup file: change the parameters as seen fit, then run this script 
% to create the inputs.dat file that will be read by emp2d. 

function [in,jobid] = emp1Drun(in)

loadconstants;

% number of time steps to run
in.tsteps = floor(3 * in.maxalt / vp / in.dt);

% need x vector; this won't be save in the output file, but we need it to
% create an ionosphere.
in.xx = in.stepalt/in.dx1 + (in.maxalt - in.stepalt)/in.dx2 + 1;
x = zeros(in.xx,1);
dx = zeros(in.xx-1,1);
x(1) = 0;
for i = 2:in.xx,
    if x(i-1) < in.stepalt,
        x(i) = x(i-1) + in.dx1;
    else
        x(i) = x(i-1) + in.dx2;
    end
    dx(i-1) = x(i) - x(i-1);
end
        
fprintf('Grid is %d cells, and will run %d time steps\n',length(x),in.tsteps);


%% input current pulse (vs time) for lightning strike. This can be anything, as long as it has tsteps values. Here, we implement some filtering to keep dr < 10*lambda_max. 

Ein = zeros(in.tsteps,1);

% okay, set up current waveform: linear rise, exp decay
for t = 1:in.tsteps,
    if (t * in.dt < in.taur),
        Ein(t) = in.E0 * t * in.dt / in.taur;
    else
        Ein(t) = in.E0 * exp(-(t*in.dt-in.taur)^2/in.tauf^2);
    end
end

% okay, filter it

b = fir1(40,in.fcut/(1/in.dt/2));

Ein2 = filter(b,1,Ein);
Ein2 = Ein2 / max(Ein2) * in.E0;   % rescaled to peak current

in.Ein = Ein2;
in.Ein(Ein < 1e-6) = 1e-6;

in.tsteps = length(in.Ein);


% set up probe points
[i,alts] = find(in.probealt < in.maxalt);
in.probex = zeros(length(i),1);
for m = 1:length(i),
    in.probex(m) = find(x >= in.probealt(alts(m)),1,'first');
end
disp(in.probex)

% write everything to a file, to be read by program

if ~exist(in.rundir,'dir'),
    mkdir(in.rundir);
end

fid = fopen([in.rundir '/inputs.dat'],'w');
fwrite(fid,in.doionosphere,'int');
fwrite(fid,in.doioniz,'int');
fwrite(fid,in.maxalt,'double');
fwrite(fid,in.stepalt,'double');
fwrite(fid,in.dx1,'double');
fwrite(fid,in.dx2,'double');
fwrite(fid,in.dt,'double');
fwrite(fid,in.tsteps,'int');
fwrite(fid,in.Ein,'double');
fwrite(fid,in.numfiles,'int');
fwrite(fid,length(in.probex),'int');
fwrite(fid,in.probex,'int');
fclose(fid);


%% magnetic field specified over 2D domain.

% vectors for magnetic field; they will be hh size, then will be made into
% matrices inside code.
in.Bx = in.Bmag(1); % cos(linspace(0,pi/2,hh));
in.By = in.Bmag(2); %sin(linspace(0,pi/2,hh));
in.Bz = in.Bmag(3); %sin(linspace(0,pi/2,hh));

fid = fopen([in.rundir '/B0.dat'],'w');
fwrite(fid,in.Bx,'double');
fwrite(fid,in.By,'double');
fwrite(fid,in.Bz,'double');
fclose(fid);


%% ionosphere and atmosphere densities. Run setupAtmosphere and save ne.dat and nd.dat.

beta = 0.8;
hk = 82;

ne = IRIionosphere1(x/1000);
nd = MSISatmosphere1(x/1000);
ndt = nd.total * 1e6;

wpmax = sqrt(max(ne)*QE^2/ME/e0);
dtmax = 1/sqrt((vp/in.dx2)^2 + (wpmax/4)^2);

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

rates = getNonlinearRates(x/1000,nd);

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


%% quick plots

h1 = figure(1);
set(h1,'position',[100 100 600 600]);
ax1 = subplot(221);
plot(ax1,log10(ne),x/1000);
hold(ax1,'on');
plot(ax1,log10(ni),x/1000,'r');
legend(ax1,'electron density','+ion density');


mue = 1.4856 * ndt(1) ./ ndt;
nue = (QE / ME) ./ mue;
nui = nue / 100;

ax2 = subplot(222);
plot(ax2,log10(nue),x/1000);
hold(ax2,'on');
plot(ax2,log10(nui),x/1000,'r');
legend(ax2,'electron coll. freq.','ion coll. freq.');

% source

ax3 = subplot(223);
tvec = 0:in.dt:(length(in.Ein)-1)*in.dt;
plot(ax3,tvec*1e6,in.Ein/1000);
set(ax3,'xlim',[0 in.tsteps*2*in.dt*1e6]);
xlabel(ax3,'time (us)');
ylabel(ax3,'Current (kA)');

ax4 = subplot(224);
plot(ax4,nd.temp,x/1000);
xlabel(ax4,'Ambient Temperature');

drawnow;


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
