% ricker wavelet source and circular polarization

clear all; close all;

dt = 1e-10;
t = 1:1000;

f0 = 150e6;

TD = sqrt(6)/pi/f0;
TR = TD/sqrt(3);
offset = 1.6*TD;

ricker = (1 - 2*(pi*f0*(t*dt-offset)).^2) .* exp(-(pi*f0*(t*dt-offset)).^2);

% create 90-degree out-of-phase ricker

offset2 = 1.6*TD + TR/2;
ricker2 = (1 - 2*(pi*f0*(t*dt-offset2)).^2) .* exp(-(pi*f0*(t*dt-offset2)).^2);

h1 = figure(1);

subplot(121);
plot(t*dt*1e9,ricker);
hold on;
plot(t*dt*1e9,ricker2,'r');

%% plot frequency response; how far above f0 do we go?

rickf = fftshift(fft(ricker));

f = linspace(-1/2/dt,1/2/dt,length(t));

subplot(122);
plot(f/1e6,abs(rickf));

set(gca,'xlim',[-500 500]);