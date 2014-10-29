%% plot 2D elve movie

clear all; close all;

datadir = '/shared/users/ram80/empcodes/runs/cidtest/I0_2e05/';

skiptotheend = 0;

domovie = 1;

datatype = 'double';

vp = 3e8;

s = get2drunparams(datadir,datatype);

eti = (s.camdist - s.range)/vp;
etf = (s.camdist + s.range)/vp + s.tsteps*s.dt;
elvedt = (etf - eti)/(s.elvesteps-1);

h1 = figure(1);
set(h1,'position',[100 100 1500 800]);
set(h1,'PaperSize',[11 8.5],'PaperPosition',[0.1 0.1 10.8 8.3]);
set(h1,'color',[1 1 1]);
set(h1,'Renderer','zbuffer');

ax1 = subplot(221);
ax2 = subplot(222);
ax3 = subplot(223);
ax4 = subplot(224);
%ax5 = subplot(325);

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
    
if domovie,
    aviobj = VideoWriter([datadir 'elvemovie.avi']);
    aviobj.FrameRate = 20;
    open(aviobj);
end

for t = 300:500,
    
    imagesc(azvec,elvec,elveN21P(:,:,t)'/1e6,'parent',ax1); axis(ax1,'xy'); caxis(ax1, [0 2]);
    imagesc(azvec,elvec,elveN22P(:,:,t)'/1e6,'parent',ax2); axis(ax2,'xy'); caxis(ax2, [0 2]);
    imagesc(azvec,elvec,elveN2P1N(:,:,t)'/1e6,'parent',ax3); axis(ax3,'xy'); caxis(ax3, [0 0.1]);
    imagesc(azvec,elvec,elveN2PM(:,:,t)'/1e6,'parent',ax4); axis(ax4,'xy'); caxis(ax4, [0 0.1]);
    %imagesc(azvec,elvec,elveO2P1N(:,:,t)'/1e6,'parent',ax5); axis(ax5,'xy'); caxis(ax5, [0 1e5]);
    title(ax1,sprintf('N2 1P at %.2f ms',t*elvedt*1e3));
    title(ax2,sprintf('N2 2P at %.2f ms',t*elvedt*1e3));
    title(ax3,sprintf('N2+ 1N at %.2f ms',t*elvedt*1e3));
    title(ax4,sprintf('N2+ M at %.2f ms',t*elvedt*1e3));
    colorbar('peer',ax1);
    colorbar('peer',ax2);
    colorbar('peer',ax3);
    colorbar('peer',ax4);
    
    drawnow;
    
    if domovie,
        F = getframe(h1);
        writeVideo(aviobj,F);
    end
    
end

if domovie,
    close(aviobj);
end


%%

imagesc(azvec,elvec,mean(elveN21P,3)','parent',ax1);
axis(ax1,'xy');
colorbar('peer',ax1);
imagesc(azvec,elvec,mean(elveN22P,3)','parent',ax2);
axis(ax2,'xy');
colorbar('peer',ax2);
imagesc(azvec,elvec,mean(elveN2P1N,3)','parent',ax3);
axis(ax3,'xy');
colorbar('peer',ax3);
imagesc(azvec,elvec,mean(elveN2PM,3)','parent',ax4);
axis(ax4,'xy');
colorbar('peer',ax4);
imagesc(azvec,elvec,mean(elveO2P1N,3)','parent',ax5);
axis(ax5,'xy');
colorbar('peer',ax5);

%% plot PIPER view


PVelveRed = zeros(16,s.elvesteps);
PHelveRed = zeros(s.elvesteps,16);
PVelveBlue = zeros(16,s.elvesteps);
PHelveBlue = zeros(s.elvesteps,16);

elvet = 0:elvedt:(s.elvesteps-1)*elvedt;

for t = 1:s.elvesteps,
    for m = 1:16,
        
        PVelveRed(m,t) = mean(mean(elveN21P(:,(m-1)*4+1:m*4,t))+mean(elveN2PM(:,(m-1)*4+1:m*4,t)));
        PVelveBlue(m,t) = mean(mean(elveN22P(:,(m-1)*4+1:m*4,t))+mean(elveN2P1N(:,(m-1)*4+1:m*4,t)));
        
    end
    
    for m = 1:32,
        PHelveRed(t,m) = mean(mean(elveN21P((m-1)*4+1:m*4,:,t))+mean(elveN2PM((m-1)*4+1:m*4,:,t)));
        PHelveBlue(t,m) = mean(mean(elveN22P((m-1)*4+1:m*4,:,t))+mean(elveN2P1N((m-1)*4+1:m*4,:,t)));
        
    end
end

h2 = figure(2);
set(h2,'position',[2300 150 800 700]);

ax21 = subplot(221);
ax22 = subplot(222);
ax23 = subplot(223);
ax24 = subplot(224);

imagesc(elvet*1e3,1:16,PVelveRed,'parent',ax21); 
axis(ax21,'xy');
xlabel(ax21,'Time (ms)');
ylabel(ax21,'Channel number');
title(ax21,'N2 1P Vertical array');

imagesc(-7:24,elvet*1e3,PHelveRed,'parent',ax22);
ylabel(ax22,'Time (ms)');
xlabel(ax22,'Channel number');
title(ax22,'N2 1P Horizontal array');
set(ax22,'xlim',[1 16]);

imagesc(elvet*1e3,1:16,PVelveBlue,'parent',ax23); axis(ax23,'xy');
xlabel(ax23,'Time (ms)');
ylabel(ax23,'Channel number');
title(ax23,'N2 2P Vertical array');

imagesc(-7:24,elvet*1e3,PHelveBlue,'parent',ax24);
ylabel(ax24,'Time (ms)');
xlabel(ax24,'Channel number');
title(ax24,'N2 2P Horizontal array');
set(ax24,'xlim',[1 16]);
