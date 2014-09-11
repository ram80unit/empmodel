% plot Probe signals from 2D run

clear all; close all;

datadir = '/shared/users/ram80/empcodes/runs/cidtest/I0_1e05/';

plotconductivity = 0;

datatype = 'double';

re = 6370e3;
qe = 1.6e-19;
me = 9.11e-31;
u0 = 4*pi*1e-7;
B0 = 10000e-9;      % fix this to read from file!
wce0 = qe*B0/me;

s = get2drunparams(datadir,datatype);

fid = fopen([datadir 'sferic.dat'],'r');
tsteps = fread(fid,1,'int');
rr = fread(fid,1,'int');
hh = fread(fid,1,'int');
numfiles = fread(fid,1,'int');
dt = fread(fid,1,datatype);
r = fread(fid,rr,datatype);
th = fread(fid,hh,datatype);
fclose(fid);

t = 0:dt:(tsteps-1)*dt;
numf = 2^(nextpow2(tsteps)+1);
f = linspace(-1/dt/2,1/dt/2,numf);

fid = fopen([datadir 'Probe.dat'],'r');
nprobes = fread(fid,1,'int');
prober = fread(fid,nprobes,'int');
probet = fread(fid,nprobes,'int');
Erprobe = fread(fid,[nprobes tsteps],datatype);
Etprobe = fread(fid,[nprobes tsteps],datatype);
Epprobe = fread(fid,[nprobes tsteps],datatype);
Hrprobe = fread(fid,[nprobes tsteps],datatype);
Htprobe = fread(fid,[nprobes tsteps],datatype);
Hpprobe = fread(fid,[nprobes tsteps],datatype);
fclose(fid);

h1 = figure(1);
set(h1,'position',[500 200 1200 700]);

ax1 = subplot(221);
ax2 = subplot(222);
ax3 = subplot(223);
ax4 = subplot(224);

plot(ax1,t*1000,Erprobe(:,:),'b'); hold(ax1,'on');
plot(ax1,t*1000,Etprobe(:,:),'r');
plot(ax1,t*1000,Epprobe(:,:),'k');
%set(ax1,'xlim',[1.65 1.9]);
xlabel(ax1,'Time (ms)');
ylabel(ax1,'V/m');
title(ax1,'E field components');

plot(ax2,t*1000,Hrprobe(:,:)*u0*1e9,'b'); hold(ax2,'on');
plot(ax2,t*1000,Htprobe(:,:)*u0*1e9,'r');
plot(ax2,t*1000,Hpprobe(:,:)*u0*1e9,'k');
%set(ax2,'xlim',[1.65 1.9]);
xlabel(ax2,'Time (ms)');
ylabel(ax2,'nT');
title(ax2,'B field components');

plot(ax3,f/1000,abs(fftshift(fft(Erprobe(:,:)',numf))),'b'); hold(ax3,'on');
plot(ax3,f/1000,abs(fftshift(fft(Etprobe(:,:)',numf))),'r');
plot(ax3,f/1000,abs(fftshift(fft(Epprobe(:,:)',numf))),'k');
set(ax3,'xlim',[0 100]);
xlabel(ax3,'Frequency (kHz)');
ylabel(ax3,'Arbitrary');
title(ax3,'E field spectra');
legend(ax3,'Er','Et','Ep');

plot(ax4,f/1000,abs(fftshift(fft(Hrprobe(:,:)',numf))),'b'); hold(ax4,'on');
plot(ax4,f/1000,abs(fftshift(fft(Htprobe(:,:)',numf))),'r');
plot(ax4,f/1000,abs(fftshift(fft(Hpprobe(:,:)',numf))),'k');
set(ax4,'xlim',[0 100]);
xlabel(ax4,'Frequency (kHz)');
ylabel(ax4,'Arbitrary');
title(ax4,'B field spectra');
legend(ax4,'Br','Bt','Bp');
