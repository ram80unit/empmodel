% mean free path from number density, temperature, etc

clear all; close all;

QE = 1.6e-19;           % electron charge
ME = 9.1e-31;           % electron mass in kg
MN = 4.65e-26;          % mass of N2 in kg
KB = 1.38e-23;          % Boltzmann constant
NA = 6.02e23;           % Avagadro's number

% altitudes of interest
h = 0:1:150;

% load atmosphere
nd = MSISatmosphere1(h);
N = nd.total * 1e6;     % total number density
T = nd.temp;            % temperature

% mean free path derived from viscosity. see Wikipedia page on MFP

vis = 18.27e-6 * (291.15 + 120) ./ (T + 120) .* (T/291.15).^1.5;
p = N*KB.*T;

lambda = vis./N .* sqrt(pi/2/MN/KB./T);


%% plot it (lazy style)

h1 = figure(1);
ax = axes;
%semilogx(ax,N,h,'k');
%hold(ax,'on');
semilogx(ax,lambda,h,'m');
%semilogx(ax,T,h,'k');

%legend('number density (m-3)','mean free path (m)','Temperature (K)');

xlabel(ax,'Mean free path (m)');
ylabel(ax,'Altitude (km)');