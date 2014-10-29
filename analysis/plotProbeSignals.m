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

probes = getProbeFields(datadir,s,datatype);

t = 0:s.dt:(s.tsteps-1)*s.dt;
numf = 2^(nextpow2(s.tsteps)+1);
f = linspace(-1/s.dt/2,1/s.dt/2,numf);

h1 = figure(1);
set(h1,'position',[500 200 1200 700]);

ax1 = subplot(221);
ax2 = subplot(222);
ax3 = subplot(223);
ax4 = subplot(224);

plot(ax1,t*1000,probes.Er(:,:),'b'); hold(ax1,'on');
plot(ax1,t*1000,probes.Et(:,:),'r');
plot(ax1,t*1000,probes.Ep(:,:),'k');
%set(ax1,'xlim',[1.65 1.9]);
xlabel(ax1,'Time (ms)');
ylabel(ax1,'V/m');
title(ax1,'E field components');

plot(ax2,t*1000,probes.Hr(:,:)*u0*1e9,'b'); hold(ax2,'on');
plot(ax2,t*1000,probes.Ht(:,:)*u0*1e9,'r');
plot(ax2,t*1000,probes.Hp(:,:)*u0*1e9,'k');
%set(ax2,'xlim',[1.65 1.9]);
xlabel(ax2,'Time (ms)');
ylabel(ax2,'nT');
title(ax2,'B field components');

plot(ax3,f/1000,abs(fftshift(fft(probes.Er(:,:)',numf))),'b'); hold(ax3,'on');
plot(ax3,f/1000,abs(fftshift(fft(probes.Et(:,:)',numf))),'r');
plot(ax3,f/1000,abs(fftshift(fft(probes.Ep(:,:)',numf))),'k');
set(ax3,'xlim',[0 100]);
xlabel(ax3,'Frequency (kHz)');
ylabel(ax3,'Arbitrary');
title(ax3,'E field spectra');
legend(ax3,'Er','Et','Ep');

plot(ax4,f/1000,abs(fftshift(fft(probes.Hr(:,:)',numf))),'b'); hold(ax4,'on');
plot(ax4,f/1000,abs(fftshift(fft(probes.Ht(:,:)',numf))),'r');
plot(ax4,f/1000,abs(fftshift(fft(probes.Hp(:,:)',numf))),'k');
set(ax4,'xlim',[0 100]);
xlabel(ax4,'Frequency (kHz)');
ylabel(ax4,'Arbitrary');
title(ax4,'B field spectra');
legend(ax4,'Br','Bt','Bp');
