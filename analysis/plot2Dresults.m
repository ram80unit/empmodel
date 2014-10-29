% quickly read results from my 2D code

clear all; close all;

runname = 'I0_2e05';
datadir = ['/shared/users/ram80/empcodes/runs/cidtest/' runname '/'];

set(0,'DefaultFigureRenderer','zbuffer');

datatype = 'double';

domovie = 1;

loadconstants;
B0 = 50000e-9;      % fix to read from file!
wce0 = QE*B0/ME;

s = get2drunparams(datadir,datatype);

drr = floor(s.rr/s.decfactor);
dhh = floor(s.hh/s.decfactor);
if s.decfactor > 1,
    decr = s.r(2:s.decfactor:end);
    decth = s.th(2:s.decfactor:end);
else
    decr = s.r;
    decth = s.th;
end

dr = diff(decr);
dth = decth(2)-decth(1);

h1 = figure(1); set(h1,'position',[200 30 2000 900]);
set(h1,'PaperSize',[11 8.5],'PaperPosition',[0.1 0.1 10.8 8.3]);
set(h1,'color',[1 1 1]);
set(h1,'Renderer','zbuffer');

ax1 = subplot(321);
ax2 = subplot(323);
ax3 = subplot(325);
ax4 = subplot(322);
ax5 = subplot(324);
ax6 = subplot(326);

% need to setup Ek so I can compare Eeff/Ek

nd = MSISatmosphere1((s.r-s.RE)/1000);
ndt = nd.total * 1e6;
nd2d = repmat(ndt',dhh,1);


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

Er = zeros(dhh,drr);
Ez = zeros(dhh,drr);
dumnu = zeros(dhh,drr);
dumTe = zeros(dhh,drr);
S = zeros(dhh,drr);
Hp = zeros(dhh,drr);
m = 0;
t = 0:(s.tsteps-1);

range = decth*decr(1)/1000;

range2 = [-flipud(range); range]';
field2 = [-flipud(dumnu); dumnu]';

% initialize plots

im1 = pcolor(ax1,range2,(decr-s.RE)/1000,field2); shading(ax1,'flat');
im2 = pcolor(ax2,range2,(decr-s.RE)/1000,field2); shading(ax2,'flat');
im3 = pcolor(ax3,range2,(decr-s.RE)/1000,field2); shading(ax3,'flat');
im4 = pcolor(ax4,range2,(decr-s.RE)/1000,field2); shading(ax4,'flat');
im5 = pcolor(ax5,range2,(decr-s.RE)/1000,field2); shading(ax5,'flat');
im6 = pcolor(ax6,range2,(decr-s.RE)/1000,field2); shading(ax6,'flat');

colormap(ax1,c2);

caxis(ax1,[-4 6]);
caxis(ax2,[-6 4]);
caxis(ax3,[4 12]);
caxis(ax4,[-5 5]);
%caxis(ax5,[0 1e7]);
caxis(ax6,[-15 -7]);

%set(ax4,'ylim',[60 110]);
%set(ax5,'ylim',[60 110]);
%set(ax6,'ylim',[60 110]);

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
ylabel(cax3,'Log10 collisions/sec');
cax4 = colorbar('peer',ax4);
ylabel(cax4,'% Change');
cax5 = colorbar('peer',ax5);
ylabel(cax5,'number density (per m3)');
cax6 = colorbar('peer',ax6);
ylabel(cax6,'J/m3');

set(ax1,'xlim',[-250 250]);
set(ax2,'xlim',[-250 250]);

% set up vector of times of writes to file

writet = find(mod(t,floor(s.tsteps/s.numfiles)) == 0);

% movie file

if domovie, 
    aviobj = VideoWriter([datadir 'fieldsmovie.avi']);
    aviobj.FrameRate = 5;
    open(aviobj);
end

%% okay, time loop

for m = 1:30,%s.numfiles,
    
    % E file
    Er = fread(Efile,[dhh drr],datatype);
    if size(Er,1) < dhh,
        break;
    end
    Et = fread(Efile,[dhh drr],datatype);
    Ep = fread(Efile,[dhh drr],datatype);
    % Jfile
    Jer = fread(Jfile,[dhh drr],datatype);
    Jet = fread(Jfile,[dhh drr],datatype);
    Jep = fread(Jfile,[dhh drr],datatype);
    Jir = fread(Jfile,[dhh drr],datatype);
    Jit = fread(Jfile,[dhh drr],datatype);
    Jip = fread(Jfile,[dhh drr],datatype);
    Jr = Jer + Jir;
    Jt = Jet + Jit;
    Jp = Jep + Jip;
    % H file
    Hr = fread(Hfile,[dhh drr],datatype);
    Ht = fread(Hfile,[dhh drr],datatype);
    Hp = fread(Hfile,[dhh drr],datatype);
    % K file
    Eeff = fread(Kfile,[dhh drr],datatype);
    Ek = fread(Kfile,[dhh drr],datatype);
    heat = fread(Kfile,[dhh drr],datatype);
    S = fread(Kfile,[dhh drr],datatype);
    Te = fread(Kfile,[dhh drr],datatype);
    
    fprintf('Max Eeff/Ek: %.4f\n',max(max(Eeff./Ek)));
    
    % D file
    ne = fread(Dfile,[dhh drr],datatype);
    nOm = fread(Dfile,[dhh drr],datatype);
    nue = fread(Dfile,[dhh drr],datatype);
    
    % O file
    nN21P = fread(Ofile,[dhh drr],datatype);
    nN22P = fread(Ofile,[dhh drr],datatype);
    nN2p1N = fread(Ofile,[dhh drr],datatype);
    nN2pM = fread(Ofile,[dhh drr],datatype);
    nO2p1N = fread(Ofile,[dhh drr],datatype);
    nOred = fread(Ofile,[dhh drr],datatype);
    nOgrn = fread(Ofile,[dhh drr],datatype);
    

    
    % create background Te and nue from first file, far end
    
    if m == 1,
        Te0 = repmat(Te(end,:),dhh,1);
        nue0 = repmat(nue(end,:),dhh,1);
        ne0 = repmat(ne(end,:),dhh,1);
    end
    
    deltadens = 100 * (ne - ne0)./ne0;
    deltadens(isinf(deltadens)) = 0;
    
    Er2 = [flipud(Er); Er]';
    Hp2 = [flipud(Hp); Hp]';
    S2 = [flipud(S); S]';
    nue2 = [flipud(nue); nue]';
    dens2 = [flipud(deltadens); deltadens]';
    opt2 = [flipud(nN21P); nN21P]';
    heat2 = [flipud(heat); heat]';
    
    set(im1,'CData',log10(abs(Er2)));
    set(im2,'CData',log10(abs(Hp2)));
    set(im3,'CData',log10(abs(nue2)));
    set(im4,'CData',dens2);
    set(im5,'CData',opt2);
    set(im6,'CData',log10(abs(heat2)));
    
    set(ti1,'String',sprintf('Er at time %.0f us',writet(m)/10));
    set(ti2,'String',sprintf('Hp at time %.0f us',writet(m)/10));
    set(ti3,'String',sprintf('Coll. Frequency at time %.0f us',writet(m)/10));
    set(ti4,'String',sprintf('Electon Density change at time %.0f us',writet(m)/10));
    set(ti5,'String',sprintf('Optical emissions (N2 1P) at time %.0f us',writet(m)/10));
    set(ti6,'String',sprintf('Integrating heating at time %.0f us',writet(m)/10));
    
    drawnow;
    
    %print(h1,'-depsc',sprintf('fieldmovie_%02d',m));
    
    if domovie,
        F = getframe(h1);
        writeVideo(aviobj,F);
    end
    
end

if domovie,
    close(aviobj);
end


fclose('all');

xinds = round(0.9*dhh):round(0.95*dhh);

Bmax = max(Hp(xinds,11)) * 1e12 * u0;
fprintf('B-phi amplitude at Midway is %.1f pT\n',Bmax);

