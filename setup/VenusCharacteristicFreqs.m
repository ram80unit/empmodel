% test Venus Daytime ionosphere

clear all; close all;

alt = 0:1:500;

QE = 1.6e-19;
ME = 9.1e-31;
B0 = 100e-9;        % weak background magnetic field
e0 = 8.854e-12;

ne(:,1) = VenusIonosphereDay(alt,1);
ne(:,2) = VenusIonosphereDay(alt,2);
%ne(:,3) = VenusIonosphere1(alt,3);

h1 = figure(1);
set(h1,'position',[1800 100 1400 700]);
ax1 = subplot(121);
ax2 = subplot(122);

semilogx(ax1,ne,alt);
set(ax1,'xlim',[1e6 1e11],'ylim',[100 400]);

wce = QE*B0/ME;

wpe = sqrt(QE^2*ne/ME/e0);

w01 = 0.5*(-wce + sqrt(wce.^2 + wpe.^2));
w02 = w01 + wce;

semilogx(ax2,w01/2/pi/1000,alt);
hold(ax2,'on');
semilogx(ax2,w02/2/pi/1000,alt,':');
plot(ax2,[wce wce]/2/pi/1000,[alt(1) alt(end)],'k:');
set(ax2,'ylim',[100 400]);
xlabel(ax2,'Frequency (kHz)');