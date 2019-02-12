% Script to read source current .mat file
%
% Author: Caitano L. da Silva (Penn State University)
%
% Created: Jan/30/2015
% Last modification: Jan/30/2015
%
clc; clear all


%% Load .mat file
load('clds_source_current.mat')

break;

% %% Interpolate as needed                      % Uncomment this block to
% zi = 6500:10:8500;                            % interpolate current (Ic)
% ti = 0:0.5e-6:70e-6;                          % to the the resolution
% [t,z,Ic] = interpolate_current(t,z,Ic,zi,ti); % defined by zi and ti


%% Plot
fh = figure(77); clf
FS = 18;
LFS = 28;
cmap = jet;
tunit = 1e6; % i.e., us
zunit = 1; % i.e., m
Icunit = 1e-3; % i.e., kA

tplot = t*tunit;
tlimvec = [min(tplot) max(tplot)];
zplot = z*zunit;
zlimvec = [min(zplot) max(zplot)]*zunit;
Icplot = Ic*Icunit;
Iclimvec = [min(min(Icplot)) max(max(Icplot))];

h1 = imagesc(tplot,zplot,Icplot,Iclimvec);
shading interp
axis xy
colorbar
colormap(cmap)
xlabel('t (us)','Fontsize',LFS)
ylabel('z (m)','Fontsize',LFS)
title('Source Current','Fontsize',LFS)
hcb=colorbar;
ylabel(hcb,'I (kA)','Fontsize',LFS);
colormap(cmap)
xlim(tlimvec)
ylim(zlimvec)
set(gca,'TickDir','out','FontSize',FS)