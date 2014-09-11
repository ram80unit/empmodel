%% plot 2D elve movie

clear all; close all;

datadir = '/shared/users/ram80/empcodes/runs/newtest/I0_1.00e05/';

doprint = 0;

datatype = 'double';

vp = 3e8;

s = get2drunparams(datadir,datatype);

elvesteps = 1000;
eti = (s.camdist - s.range)/vp;
etf = (s.camdist + s.range)/vp + s.tsteps*s.dt;
elvedt = (etf - eti)/(elvesteps-1);

h1 = figure(1);
set(h1,'position',[100 100 1800 1000]);
set(h1,'PaperPositionMode','auto');

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

if camtype == 1,
    error('Using regular camera pixel array? Use plot2Delve.m instead');
end

fid = fopen([datadir 'elve.dat'],'r');

elveN21P = fread(fid,totalpixels*elvesteps,datatype) / elvedt;
elveN22P = fread(fid,totalpixels*elvesteps,datatype) / elvedt;
elveN2P1N = fread(fid,totalpixels*elvesteps,datatype) / elvedt;
elveN2PM = fread(fid,totalpixels*elvesteps,datatype) / elvedt;
elveO2P1N = fread(fid,totalpixels*elvesteps,datatype) / elvedt;

elveN21P(isnan(elveN21P)) = 0;
elveN22P(isnan(elveN22P)) = 0;
elveN2P1N(isnan(elveN2P1N)) = 0;
elveN2PM(isnan(elveN2PM)) = 0;
elveO2P1N(isnan(elveO2P1N)) = 0;

fclose(fid);
    
[azmat,elmat] = meshgrid(az,el);


elveN21P = reshape(elveN21P,elvesteps,totalpixels)';
elveN22P = reshape(elveN22P,elvesteps,totalpixels)';
elveN2P1N = reshape(elveN2P1N,elvesteps,totalpixels)';
elveN2PM = reshape(elveN2PM,elvesteps,totalpixels)';
elveO2P1N = reshape(elveO2P1N,elvesteps,totalpixels)';

p1 = scatter(ax1,az*180/pi,el*180/pi,150,elveN21P(:,1)/1e6,'filled');
p2 = scatter(ax2,az*180/pi,el*180/pi,150,elveN22P(:,1)/1e6,'filled');
p3 = scatter(ax3,az*180/pi,el*180/pi,150,elveN2P1N(:,1)/1e6,'filled');
p4 = scatter(ax4,az*180/pi,el*180/pi,150,elveN2PM(:,1)/1e6,'filled');
p5 = scatter(ax5,az*180/pi,el*180/pi,150,elveO2P1N(:,1)/1e6,'filled');

set([ax1 ax2 ax3 ax4 ax5],'Color',[0.5 0.5 1]);

caxis(ax1,[0 0.1*max(max(elveN21P))/1e6]);
caxis(ax2,[0 0.1*max(max(elveN22P))/1e6]);
caxis(ax3,[0 0.1*max(max(elveN2P1N))/1e6]);
caxis(ax4,[0 0.1*max(max(elveN2PM))]/1e6);
caxis(ax5,[0 0.1*max(max(elveO2P1N))]/1e6);

colorbar('peer',ax1);
colorbar('peer',ax2);
colorbar('peer',ax3);
colorbar('peer',ax4);
colorbar('peer',ax5);

ti1 = title(ax1,'N2 1P');
title(ax2,'N2 2P');
title(ax3,'N2+ 1N');
title(ax4,'N2+ M');
title(ax5,'O2+ 1N');
xlabel(ax3,'Degrees Azimuth');
ylabel(ax1,'Degrees Elevation');

for t = 300:600,
    
    set(p1,'CData',elveN21P(:,t)/1e6);
    set(p2,'CData',elveN22P(:,t)/1e6);
    set(p3,'CData',elveN2P1N(:,t)/1e6);
    set(p4,'CData',elveN2PM(:,t)/1e6);
    set(p5,'CData',elveO2P1N(:,t)/1e6);
    set(ti1,'string',sprintf('N2 1P at %.3f ms',t*elvedt*1e3));
    
    if doprint,
        print(h1,'-dpng',sprintf('AugerElve%03d.png',t));
    end
    
    drawnow;
    
end

break;

%%

imagesc(azvec,elvec,mean(elveN21P,2)','parent',ax1);
axis(ax1,'xy');
colorbar('peer',ax1);
imagesc(azvec,elvec,mean(elveN22P,2)','parent',ax2);
axis(ax2,'xy');
colorbar('peer',ax2);
imagesc(azvec,elvec,mean(elveN2P1N,3)','parent',ax3);
axis(ax3,'xy');
colorbar('peer',ax3);
imagesc(azvec,elvec,mean(elveN2PM,2)','parent',ax4);
axis(ax4,'xy');
colorbar('peer',ax4);
imagesc(azvec,elvec,mean(elveO2P1N,2)','parent',ax5);
axis(ax5,'xy');
colorbar('peer',ax5);


