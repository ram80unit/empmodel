% quickly read results from my 3D code

clear all; close all;

%runname = 'I0_100';
%datadir = ['/shared/users/ram80/empcodes/runs/heating/trans3/full/' runname '/'];
runname = 'I0_1e05';
datadir = ['/shared/users/ram80/empcodes/runs/sourcetest3/' runname '/'];

datatype = 'double';

loadconstants;
B0 = 50000e-9;      % fix to read from file!
wce0 = QE*B0/ME;

s = get3drunparams(datadir);

drr = floor(s.rr/s.decfactor);
dhh = floor(s.hh/s.decfactor);
dpp = floor(s.pp/s.decfactor);
decr = decimate(s.r,s.decfactor);
decth = decimate(s.th,s.decfactor);
decph = decimate(s.ph,s.decfactor);

dr = diff(decr);
dth = decth(2)-decth(1);
dph = decph(2)-decph(1);


h1 = figure(1); set(h1,'position',[300 30 1600 900]);
set(h1,'PaperSize',[11 7],'PaperPosition',[0.1 0.1 10.8 6.8]);
set(h1,'color',[1 1 1]);

set(h1,'renderer','painters');

ax1 = subplot(321);
ax2 = subplot(322);
ax3 = subplot(323);
ax4 = subplot(324);
ax5 = subplot(325);
ax6 = subplot(326);

nd = MSISatmosphere1((s.r-RE)/1000);
ndt = nd.total * 1e6;
nd2d = repmat(ndt',dhh,1);

fid = fopen([datadir 'ne.dat']);
ne0 = fread(fid,s.rr,datatype);
fclose(fid);
ne02d = repmat(ne0',s.hh,1);


c = colormap('jet');
c3 = [0 0 0; 0 0 0.1; 0 0 0.2; 0 0 0.3; 0 0 0.4; c];
c2 = [1 1 1; 0.8 0.8 0.9; 0.6 0.6 0.8; 0.4 0.4 0.7; 0.2 0.2 0.6; c];

% open files for reading

Efile = fopen([datadir 'output_E.dat'],'r');
Jfile = fopen([datadir 'output_J.dat'],'r');
Hfile = fopen([datadir 'output_H.dat'],'r');
Kfile = fopen([datadir 'output_K.dat'],'r');
Dfile = fopen([datadir 'output_D.dat'],'r');
Ofile = fopen([datadir 'output_O.dat'],'r');

rr = s.rr;
hh = s.hh;
pp = s.pp;

Eslicep = zeros(rr,hh);
Eslicet = zeros(rr,pp);
Eslicer = zeros(hh,pp);
Jslicep = zeros(rr,hh);
Jslicet = zeros(rr,pp);
Jslicer = zeros(hh,pp);
Hslicep = zeros(rr,hh);
Hslicet = zeros(rr,pp);
Hslicer = zeros(hh,pp);

Eeffslice = zeros(rr,hh);
Ekslice = zeros(rr,hh);
heatslice = zeros(rr,hh);
Sslice = zeros(rr,hh);
Teslice = zeros(rr,hh);
neslice = zeros(rr,hh);
n0mslice = zeros(rr,hh);
nuslice = zeros(rr,hh);
nN21Pslice = zeros(rr,hh);
nN22Pslice = zeros(rr,hh);
nN2P1Nslice = zeros(rr,hh);
nN2PMslice = zeros(rr,hh);
nO2P1Nslice = zeros(rr,hh);

t = 0:(s.tsteps-1);

range = (s.th-pi/2)*s.r(s.nground+1)/1000;

newr = min(s.r):1000:max(s.r);

newE = interp2(range,s.r',Eslicep,range,newr,'linear');
newS = interp2(range,s.r',Sslice,range,newr,'linear');
newnu = interp2(range,s.r',nuslice,range,newr,'linear');
newne = interp2(range,s.r',neslice,range,newr,'linear');
newJ = interp2(range,s.r',Jslicep,range,newr,'linear');
newTe = interp2(range,s.r',Teslice,range,newr,'linear');

% initialize plots

im1 = imagesc(range,(newr-s.RE)/1000,log10(abs(newE)+1e-8),'parent',ax1); axis(ax1,'xy');
im2 = imagesc(range,(newr-s.RE)/1000,log10(abs(newS)),'parent',ax2); axis(ax2,'xy');
im3 = imagesc(range,(newr-s.RE)/1000,newnu,'parent',ax3); axis(ax3,'xy');
im4 = imagesc(range,(newr-s.RE)/1000,newne,'parent',ax4); axis(ax4,'xy');
im5 = imagesc(range,(newr-s.RE)/1000,log10(abs(newJ)),'parent',ax5); axis(ax5,'xy');
im6 = imagesc(range,(newr-s.RE)/1000,newTe,'parent',ax6); axis(ax6,'xy');

colormap(ax1,c2);

caxis(ax1,[-5 3]);
caxis(ax2,[-10 -2]);
%caxis(ax3,[0 5000]);
caxis(ax4,[-1 1]);
caxis(ax5,[-12 -4]);
caxis(ax6,[0 100]);

xlabel(ax5,'Range along ground (km)');
ylabel(ax1,'Altitude (km)');
ylabel(ax2,'Altitude (km)');
ylabel(ax3,'Altitude (km)');
ylabel(ax4,'Altitude (km)');
ylabel(ax5,'Altitude (km)');
ylabel(ax6,'Altitude (km)');
ti1 = title(ax1,'');
ti2 = title(ax2,'');
ti3 = title(ax3,'');
ti4 = title(ax4,'');
ti5 = title(ax5,'');
ti6 = title(ax6,'');

cax1 = colorbar('peer',ax1);
ylabel(cax1,'Log10 V/m');
cax2 = colorbar('peer',ax2);
ylabel(cax2,'Log10 W/m3');
cax3 = colorbar('peer',ax3);
ylabel(cax3,'% Change');
cax4 = colorbar('peer',ax4);
ylabel(cax4,'% Change');
cax5 = colorbar('peer',ax5);
ylabel(cax5,'Log10 A/m2');
cax6 = colorbar('peer',ax6);
ylabel(cax6,'Degrees');

set(ax4,'ylim',[60 100]);
set(ax5,'ylim',[60 100]);
set(ax6,'ylim',[60 100]);


% set up vector of times of writes to file

writet = find(mod(t,floor(s.tsteps/s.numfiles)) == 0);

for m = 1:s.numfiles,
    
    Eslicep = fread(Efile,[hh rr],datatype);
    
    if size(Eslicep,1) < hh,
        fclose('all');
        break;
    end
    
    Eslicet = fread(Efile,[pp rr],datatype);
    Eslicer = fread(Efile,[pp hh],datatype);
    
    Jslicep = fread(Jfile,[hh rr],datatype);
    Jslicet = fread(Jfile,[pp rr],datatype);
    Jslicer = fread(Jfile,[pp hh],datatype);
    
    Hslicep = fread(Hfile,[hh rr],datatype);
    Hslicet = fread(Hfile,[pp rr],datatype);
    Hslicer = fread(Hfile,[pp hh],datatype);
    
    Eeffslice = fread(Kfile,[hh rr],datatype);
    Ekslice = fread(Kfile,[hh rr],datatype);
    heatslice = fread(Kfile,[hh rr],datatype);
    Sslice = fread(Kfile,[hh rr],datatype);
    Teslice = fread(Kfile,[hh rr],datatype);
    
    neslice = fread(Dfile,[hh rr],datatype);
    nOmslice = fread(Dfile,[hh rr],datatype);
    nuslice = fread(Dfile,[hh rr],datatype);

    nN21Pslice = fread(Ofile,[hh rr],datatype);
    nN22Pslice = fread(Ofile,[hh rr],datatype);
    nN2P1Nslice = fread(Ofile,[hh rr],datatype);
    nN2PMslice = fread(Ofile,[hh rr],datatype);
    nO2P1Nslice = fread(Ofile,[hh rr],datatype);
    
    %fprintf('Max Eeff/Ek: %.4f\n',max(max(Eeffslice./Ek)));
    
    deltadens = 100 * (neslice - ne02d)./ne02d;
    deltadens(isinf(deltadens)) = 0;
    
    if m == 1,
        Te0 = Teslice;
        nu0 = nuslice;
    end
    deltaTe = Teslice - Te0;
    deltanu = 100 * (nuslice - nu0) ./ nu0;
    
    newE = interp2(range,s.r',Eslicet',range,newr,'linear');
    newS = interp2(range,s.r',Sslice',range,newr,'linear');
    newnu = interp2(range,s.r',deltanu',range,newr,'linear');
    newne = interp2(range,s.r',deltadens',range,newr,'linear');
    newJ = interp2(range,s.r',Jslicep',range,newr,'linear');
    newTe = interp2(range,s.r',deltaTe',range,newr,'linear');
    
    set(im1,'CData',log10(abs(newE)));
    set(im2,'CData',log10(abs(newS)));
    set(im3,'CData',newnu);
    set(im4,'CData',newne);
    set(im5,'CData',log10(abs(newJ)));
    set(im6,'CData',newTe);
    
    set(ti1,'String',sprintf('Er at time %.0f us',writet(m)/10));
    set(ti2,'String',sprintf('Poyting Flux S at time %.0f us',writet(m)/10));
    set(ti3,'String',sprintf('%% change in nu at time %.0f us',writet(m)/10));
    set(ti4,'String',sprintf('%% change in ne at time %.0f us',writet(m)/10));
    set(ti5,'String',sprintf('Jr at time %.0f us',writet(m)/10));
    set(ti6,'String',sprintf('Change in Te (K) at time %.0f us',writet(m)/10));
    
    cmax = max([20 max(max(deltaTe))]);
    set(ax6,'clim',[0 cmax]);
    
    drawnow;
    pause(0.2);
    
end

fclose('all');


