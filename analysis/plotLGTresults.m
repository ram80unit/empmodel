% quickly read results from my 2D code

clear all; close all;

datadir = '/shared/users/ram80/empcodes/lightning/';

datatype = 'double';
printname = 'upward';

re = 6370e3;
qe = 1.6e-19;
me = 9.11e-31;
u0 = 4*pi*1e-7;

fid = fopen([datadir 'sferic.dat'],'r');
tsteps = fread(fid,1,'int');
rr = fread(fid,1,'int');
hh = fread(fid,1,'int');
numfiles = fread(fid,1,'int');
dt = fread(fid,1,datatype);
r = fread(fid,rr,datatype);
th = fread(fid,hh,datatype);
fclose(fid);

dr = diff(r);
dth = th(2)-th(1);

h1 = figure(1); set(h1,'position',[300 30 1100 900]);
set(h1,'PaperSize',[11 7],'PaperPosition',[0.1 0.1 10.8 6.8]);
set(h1,'color',[1 1 1]);

ax1 = subplot(221);
ax2 = subplot(222);
ax3 = subplot(223);

c = colormap('jet');
c2 = [1 1 1; 0.8 0.8 0.9; 0.6 0.6 0.8; 0.4 0.4 0.7; 0.2 0.2 0.6; c];

% open files for reading

Efile = fopen([datadir 'output_fields.dat'],'r');

Er = zeros(hh,rr);
Et = zeros(hh,rr);
Hp = zeros(hh,rr);

t = 0:(tsteps-1);
range = th*r(1)/1000;

% initialize plots

im1 = pcolor(range,(r-re)/1000,log10(abs(Er')+1e-8),'parent',ax1); shading(ax1,'flat');
im2 = pcolor(range,(r-re)/1000,log10(abs(Et')+1e-8),'parent',ax2); shading(ax2,'flat');
im3 = pcolor(range,(r-re)/1000,log10(abs(Hp')*377+1e-8),'parent',ax3); shading(ax3,'flat');

ti1 = title(ax1,'');
ti2 = title(ax2,'');
ti3 = title(ax3,'');

colormap(ax1,c2);

caxis(ax1,[-2 5]);
caxis(ax2,[-2 5]);
caxis(ax3,[-2 5]);

xlabel(ax1,'Range (km)');
ylabel(ax1,'Altitude (km)');
xlabel(ax2,'Range (km)');
ylabel(ax2,'Altitude (km)');
xlabel(ax3,'Range (km)');
ylabel(ax3,'Altitude (km)');

%% time loop

for m = 1:numfiles,
    
    Er = fread(Efile,[hh rr],datatype);
    
    if size(Er,1) < hh,
        break;
    end
    
    Et = fread(Efile,[hh rr],datatype);
    Hp = fread(Efile,[hh rr],datatype);
    
    set(im1,'CData',log10(abs(Er)'));
    set(im2,'CData',log10(abs(Et)'));
    set(im3,'CData',log10(abs(Hp)'*377));
    set(ti1,'String',sprintf('Er at t = %.0f us',tsteps*dt*m/100*1e6));
    set(ti2,'String',sprintf('Et at t = %.0f us',tsteps*dt*m/100*1e6));
    set(ti3,'String',sprintf('Hp at t = %.0f us',tsteps*dt*m/100*1e6));
    
    drawnow;
    
    %print(h1,'-dpdf',sprintf('%s/figs/%s%03d.pdf',datadir,printname,m));

    
end

%% plot probe fields

fid = fopen([datadir 'Probe.dat'],'r');
Erprobe = fread(fid,tsteps,datatype);
Etprobe = fread(fid,tsteps,datatype);
Hpprobe = fread(fid,tsteps,datatype);
fclose(fid);

ax4 = subplot(224);
t = (1:tsteps)*dt*1e6;
plot(ax4,t,Erprobe,'k');
hold(ax4,'on');
plot(ax4,t,Etprobe,'b');
plot(ax4,t,Hpprobe*377,'r');
legend('Er','Et','Hp * eto0');
xlabel(ax4,'Time (us)');
ylabel(ax4,'Field amplitude (V/m)');
title(ax4,'Fields measured at 37 km along ground');
set(ax4,'xlim',[80 300]);

%print(h1,'-dpdf',sprintf('%s/figs/%s%03d.pdf',datadir,printname,m+1));
