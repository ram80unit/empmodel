%% Read BOLSIG+ outputs for air from cd

function [EoverN, outputs] = readoutputs_air(bolsigfile,doplot,reldens)

% simulation:
% 100 E/N bins, from 0.04 Td to 4000 Td
% All O, He, Ar, N2 and O2 cross sections
%
% upon clicking "save data", I chose "Separate Tables vs E/N" and uncheck
% "selected items only". First three boxes on left are checked.

fid = fopen(bolsigfile,'rt');

% first section is constants for all runs

for m = 1:20,
    dummy = fgetl(fid);
end

% First set of data

for k = 1:20,
    
    fields(k).title = fgetl(fid);
    data = fscanf(fid,'%f %f\n',[2 100]);
    fields(k).var1 = data(1,:);
    fields(k).var2 = data(2,:);
    
end

% second set, has an extra title line

for k = 21:74,
    
    fields(k).title1 = fgetl(fid);
    fields(k).title2 = fgetl(fid);
    data = fscanf(fid,'%f %f\n',[2 100]);
    fields(k).var1 = data(1,:);
    fields(k).var2 = data(2,:);
    
end

fclose(fid);

%% plot a few important parameters: total ionization, attachment, mobility, collisions

% O ionization: number 25
% N2 ionization: 51
% O2 ionization: 66
% O2 attachment: 67 and 68
% Ar ionization: 71
% He ionization: 74

N0 = 2.688e25;
Tdfac = 1e-21;
EoverN = fields(1).var1 * Tdfac;
EN0overN = EoverN * N0;

outputs.ioniz = reldens.o * fields(25).var2 + reldens.n2 * fields(51).var2 + ...
    reldens.o2 * fields(66).var2 + reldens.ar * fields(71).var2 + ...
    reldens.he * fields(74).var2;

% attachment is numbers 67 and 68

outputs.attach = reldens.o2 * (fields(67).var2 + fields(68).var2);

% Mobility is number 4

outputs.mobility = fields(4).var2;

% Electron mean energy is number 3

outputs.energy = fields(3).var2;

outputs.Ored = reldens.o * fields(22).var2;
outputs.Ogrn = reldens.o * fields(23).var2;
outputs.Oblk = reldens.o * fields(24).var2;

outputs.n21p = reldens.n2 * fields(40).var2;
outputs.n22p = reldens.n2 * fields(47).var2;

if doplot,
    
    he = figure;
    set(he,'position',[300 100 900 800]);
    
    ax1 = subplot(221);
    loglog(ax1, EN0overN, outputs.ioniz * N0, 'bx-'); hold on;
    title(ax1, 'Total Ionization');
    xlabel(ax1, 'E N_0 / N (V/m)');
    ylabel(ax1, 'v_i N_0 (1/s)');
    set(ax1, 'xlim', [5e5 6e7]);
    set(ax1, 'ylim', [1e4 1e13]);
    
    ax2 = subplot(222);
    loglog(ax2, EN0overN, outputs.attach * N0, 'bx-'); hold on;
    title(ax2, 'Two-body Attachment');
    xlabel(ax2, 'E N_0 / N (V/m)');
    ylabel(ax2, 'v_a N_0 (1/s)');
    set(ax2, 'xlim', [5e5 6e7]);
    set(ax2, 'ylim', [1e4 1e9]);
    
    ax3 = subplot(223);
    loglog(ax3, EN0overN, outputs.mobility / N0, 'bx-'); hold on;
    title(ax3, 'Electron Mobility');
    xlabel(ax3, 'E N0 / N (V/m)');
    ylabel(ax3, '\mu_e N/N_0 (m^2/V/s)');
    set(ax3, 'xlim', [1e3 1e8]);
    set(ax3, 'ylim', [1e-2 1e1]);
    
    ax4 = subplot(224);
    loglog(ax4, EN0overN, outputs.energy, 'bx-'); hold on;
    title(ax4, 'Electron mean Energy');
    xlabel(ax4, 'E N0 / N (V/m)');
    ylabel(ax4, 'Mean energy (eV)');
    set(ax4, 'xlim', [1e3 1e8]);
    set(ax4, 'ylim', [1e-2 1e2]);
    
end
%% Plot Elendif outputs from Moss on top....

e = EN0overN;
alt = 0;

for m = 1:length(e),
    outputs.mossattach(m) = air1(e(m),alt,2);
    outputs.mossioniz(m) = air1(e(m),alt,10);
    outputs.mossenergy(m) = air1(e(m),alt,1);
    outputs.mossmobility(m) = air1(e(m),alt,11);
    outputs.mossn21p(m) = air1(e(m),alt,5) / N0;
    outputs.mossn22p(m) = air1(e(m),alt,6) / N0;
    outputs.mossn2p1n(m) = air1(e(m),alt,7) / N0;
    outputs.mosso2p1n(m) = air1(e(m),alt,8) / N0;
    outputs.mossn2pm(m) = air1(e(m),alt,9) / N0;
end

if doplot,
    
    loglog(ax1, e, outputs.mossioniz, 'rx-');
    loglog(ax2, e, outputs.mossattach, 'rx-');
    loglog(ax3, e, outputs.mossmobility, 'rx-');
    loglog(ax4, e, outputs.mossenergy, 'rx-');
    
end
