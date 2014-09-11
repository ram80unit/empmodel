%% plot 3D elve movie

clear all; close all;

datadir = '/shared/users/ram80/empcodes/runs/elvenogw/gwavemag_0.2/';

skiptotheend = 1;

datatype = 'double';

vp = 3e8;

s = get3drunparams(datadir);

eti = (s.camdist - s.range)/vp;
etf = (s.camdist + s.range)/vp + s.tsteps*s.dt;
elvedt = (etf - eti)/(s.elvesteps-1);

h1 = figure(1);
set(h1,'position',[100 100 1200 800]);

ax1 = subplot(321);
ax2 = subplot(322);
ax3 = subplot(323);
ax4 = subplot(324);
ax5 = subplot(325);

fid = fopen([datadir 'camera.dat'],'r');
camtype = fread(fid,1,'int');
totalpixels = fread(fid,1,'int');
az = fread(fid,totalpixels,'double');
el = fread(fid,totalpixels,'double');
fclose(fid);

if camtype == 0,
    error('Using Auger pixel array? Use plot2DelveAuger.m instead');
end

numaz = length(unique(az));
numel = length(unique(el));
imsize = totalpixels;

fid = fopen([datadir 'elve.dat'],'r');

elveN21P = fread(fid,imsize*s.elvesteps,datatype) / elvedt;
elveN22P = fread(fid,imsize*s.elvesteps,datatype) / elvedt;
elveN2P1N = fread(fid,imsize*s.elvesteps,datatype) / elvedt;
elveN2PM = fread(fid,imsize*s.elvesteps,datatype) / elvedt;
elveO2P1N = fread(fid,imsize*s.elvesteps,datatype) / elvedt;

elveN21P(isnan(elveN21P)) = 0;
elveN22P(isnan(elveN22P)) = 0;
elveN2P1N(isnan(elveN2P1N)) = 0;
elveN2PM(isnan(elveN2PM)) = 0;
elveO2P1N(isnan(elveO2P1N)) = 0;

fclose(fid);
    
elveN21P = permute(reshape(elveN21P,s.elvesteps,numel,numaz),[3 2 1]);
elveN22P = permute(reshape(elveN22P,s.elvesteps,numel,numaz),[3 2 1]);
elveN2P1N = permute(reshape(elveN2P1N,s.elvesteps,numel,numaz),[3 2 1]);
elveN2PM = permute(reshape(elveN2PM,s.elvesteps,numel,numaz),[3 2 1]);
elveO2P1N = permute(reshape(elveO2P1N,s.elvesteps,numel,numaz),[3 2 1]);

azvec = unique(az)*180/pi;
elvec = unique(el)*180/pi;
    

for t = 1:s.elvesteps,
    
    imagesc(azvec,elvec,elveN21P(:,:,t)','parent',ax1); axis(ax1,'xy'); caxis(ax1, [0 3e6]);
    imagesc(azvec,elvec,elveN22P(:,:,t)','parent',ax2); axis(ax2,'xy'); caxis(ax2, [0 1e6]);
    imagesc(azvec,elvec,elveN2P1N(:,:,t)','parent',ax3); axis(ax3,'xy'); caxis(ax3, [0 1e5]);
    imagesc(azvec,elvec,elveN2PM(:,:,t)','parent',ax4); axis(ax4,'xy'); caxis(ax4, [0 1e5]);
    imagesc(azvec,elvec,elveO2P1N(:,:,t)','parent',ax5); axis(ax5,'xy'); caxis(ax5, [0 1e5]);
    title(ax1,sprintf('%.2f ms',t*elvedt*1e3));
    
    drawnow;
    
end

%%

imagesc(azvec,elvec,mean(elveN21P,3)'/1e6,'parent',ax1);
axis(ax1,'xy');
ylabel(ax1,'Elevation (degrees)'); title(ax1,'Time-integrated elve in N2 1P bands');
colorbar('peer',ax1);

imagesc(azvec,elvec,mean(elveN22P,3)'/1e6,'parent',ax2);
axis(ax2,'xy');
title(ax2,'N2 2P bands');
colorbar('peer',ax2);

imagesc(azvec,elvec,mean(elveN2P1N,3)'/1e6,'parent',ax3);
axis(ax3,'xy');
ylabel(ax3,'Elevation (degrees)'); title(ax3,'N2+ 1N bands');
colorbar('peer',ax3);

imagesc(azvec,elvec,mean(elveN2PM,3)'/1e6,'parent',ax4);
axis(ax4,'xy');
xlabel(ax4,'Azimuth (degrees)'); title(ax4,'N2+ M bands');
colorbar('peer',ax4);

imagesc(azvec,elvec,mean(elveO2P1N,3)'/1e6,'parent',ax5);
axis(ax5,'xy');
xlabel(ax5,'Azimuth (degrees)'); title(ax5,'O2+ 1N bands');
colorbar('peer',ax5);

% color elve! Just for fun and presentations
% need to modify this to take proper RGB filters into account

elvered = mean(elveN21P,3)' + mean(elveN2PM,3)' + mean(elveO2P1N,3)';
elveblue = mean(elveN22P,3)' + mean(elveN2P1N,3)';

elvecolor = zeros(size(elvered,1),size(elvered,2),3);
elvecolor(:,:,1) = elvered/max(max(elvered));
elvecolor(:,:,3) = elveblue/max(max(elvered));

h2 = figure(2);
set(h2,'position',[700 100 1000 500]);
axh2 = axes;

imagesc(azvec,elvec,elvecolor,'parent',axh2);
axis(axh2,'xy');
xlabel(axh2,'Azimuth (degrees)');
ylabel(axh2,'Elevation (degrees)');
title(axh2,'Simulated Real-color Elve');


%% plot PIPER view


PVelveRed = zeros(16,s.elvesteps);
PHelveRed = zeros(s.elvesteps,16);
PVelveBlue = zeros(16,s.elvesteps);
PHelveBlue = zeros(s.elvesteps,16);

elvet = 0:elvedt:(s.elvesteps-1)*elvedt;

for t = 1:s.elvesteps,
    for m = 1:16,
    
    PVelveRed(m,t) = mean(mean(elveN21P(:,(m-1)*4+1:m*4,t))+mean(elveN2PM(:,(m-1)*4+1:m*4,t)));
    PHelveRed(t,m) = mean(mean(elveN21P((m-1)*8+1:m*8,:,t))+mean(elveN2PM((m-1)*8+1:m*8,:,t)));
    PVelveBlue(m,t) = mean(mean(elveN22P(:,(m-1)*4+1:m*4,t))+mean(elveN2P1N(:,(m-1)*4+1:m*4,t)));
    PHelveBlue(t,m) = mean(mean(elveN22P((m-1)*8+1:m*8,:,t))+mean(elveN2P1N((m-1)*8+1:m*8,:,t)));
    
    end
end

h3 = figure(3);
set(h3,'position',[1700 150 800 700]);

ax21 = subplot(221);
ax22 = subplot(222);
ax23 = subplot(223);
ax24 = subplot(224);

imagesc(elvet*1e3,1:16,PVelveRed,'parent',ax21); 
axis(ax21,'xy');
xlabel(ax21,'Time (ms)');
ylabel(ax21,'Channel number');
title(ax21,'PIPER: N2 1P Vertical array');

imagesc(1:16,elvet*1e3,PHelveRed,'parent',ax22);
ylabel(ax22,'Time (ms)');
xlabel(ax22,'Channel number');
title(ax22,'PIPER: N2 1P Horizontal array');

imagesc(elvet*1e3,1:16,PVelveBlue,'parent',ax23); axis(ax23,'xy');
xlabel(ax23,'Time (ms)');
ylabel(ax23,'Channel number');
title(ax23,'PIPER: N2 2P Vertical array');

imagesc(1:16,elvet*1e3,PHelveBlue,'parent',ax24);
ylabel(ax24,'Time (ms)');
xlabel(ax24,'Channel number');
title(ax24,'PIPER: N2 2P Horizontal array');
