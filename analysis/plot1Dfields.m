% quickly read results from my 2D code

clear all; close all;

runname = 'E0_20';
datadir = ['/shared/users/ram80/empcodes/runs/test1d/' runname '/'];

datatype = 'double';

loadconstants;
B0 = 50000e-9;      % fix to read from file!
wce0 = QE*B0/ME;
Z0 = sqrt(u0/e0);

s = get1drunparams(datadir,datatype);


h1 = figure(1); set(h1,'position',[200 30 1500 800]);
set(h1,'PaperSize',[11 8.5],'PaperPosition',[0.1 0.1 10.8 8.3]);
set(h1,'color',[1 1 1]);
set(h1,'Renderer','zbuffer');

ax1 = subplot(131);
ax2 = subplot(132);
ax3 = subplot(133);

% open files for reading

Efile = fopen([datadir 'output_E.dat'],'r');
Jfile = fopen([datadir 'output_J.dat'],'r');
Hfile = fopen([datadir 'output_H.dat'],'r');

Ex = zeros(s.xx,1);
Ey = zeros(s.xx,1);
Ez = zeros(s.xx,1);
Hx = zeros(s.xx,1);
Hy = zeros(s.xx,1);
Hz = zeros(s.xx,1);
Jx = zeros(s.xx,1);
Jy = zeros(s.xx,1);
Jz = zeros(s.xx,1);

p11 = plot(ax1,log10(abs(Ex)),s.x/1000,'k');
hold(ax1,'on');
p12 = plot(ax1,log10(abs(Ey)),s.x/1000,'b');
p13 = plot(ax1,log10(abs(Ez)),s.x/1000,'r');

p21 = plot(ax2,log10(abs(Hx*Z0)),s.x/1000,'k');
hold(ax2,'on');
p22 = plot(ax2,log10(abs(Hy*Z0)),s.x/1000,'b');
p23 = plot(ax2,log10(abs(Hz*Z0)),s.x/1000,'r');

p31 = plot(ax3,Jx,s.x/1000,'k');
hold(ax3,'on');
p32 = plot(ax3,Jy,s.x/1000,'b');
p33 = plot(ax3,Jz,s.x/1000,'r');

set(ax1,'xlim',[-3 2],'ylim',[0 s.maxalt/1e3]);
set(ax2,'xlim',[-3 2],'ylim',[0 s.maxalt/1e3]);
set(ax3,'xlim',[-9 -4],'ylim',[0 s.maxalt/1e3]);

for m = 1:s.numfiles,
    
    % E file
    Ex = fread(Efile,s.xx,datatype);
    if length(Ex) < s.xx,
        break;
    end
    Ey = fread(Efile,s.xx,datatype);
    Ez = fread(Efile,s.xx,datatype);
    % Hfile
    Hx = fread(Hfile,s.xx,datatype);
    Hy = fread(Hfile,s.xx,datatype);
    Hz = fread(Hfile,s.xx,datatype);
    % Jfile
    Jex = fread(Jfile,s.xx,datatype);
    Jey = fread(Jfile,s.xx,datatype);
    Jez = fread(Jfile,s.xx,datatype);
    Jix = fread(Jfile,s.xx,datatype);
    Jiy = fread(Jfile,s.xx,datatype);
    Jiz = fread(Jfile,s.xx,datatype);
    Jx = Jex + Jix;
    Jy = Jey + Jiy;
    Jz = Jez + Jiz;
    
    set(p11,'XData',log10(abs(Ex)));
    set(p12,'XData',log10(abs(Ey)));
    set(p13,'XData',log10(abs(Ez)));
    set(p21,'XData',log10(abs(Hx*Z0)));
    set(p22,'XData',log10(abs(Hy*Z0)));
    set(p23,'XData',log10(abs(Hz*Z0)));
    set(p31,'XData',log10(abs(Jx)));
    set(p32,'XData',log10(abs(Jy)));
    set(p33,'XData',log10(abs(Jz)));

    drawnow;
    
    %print(h1,'-depsc',sprintf('fieldmovie_%02d',m));
    
end

fclose('all');

%
%Bmax = max(Hp(xinds,11)) * 1e12 * u0;
%fprintf('B-phi amplitude at Midway is %.1f pT\n',Bmax);

