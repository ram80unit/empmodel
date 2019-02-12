% 2D nd test

clear all; close all;

Re = 6370e3;

r = Re:1000:(Re+110e3);
range = -200e3:2000:200e3;
th = range/Re;
ph = range/Re;
dth = th(2) - th(1);
dph = dth;
H = 6e3;
hmax = 90e3;
gwmag = 0.5;
gwkh = 2*pi/60e3;
gwkp = 2*pi/500e3;
gwkr = 2*pi/20e3;

nd = MSISatmosphere1((r-Re)/1000);
ndt = nd.total * 1e6;

nd3 = zeros(length(r),length(th),length(ph));

for i = 1:length(r),
    for j = 1:length(th),
        for k = 1:length(ph),
            
            gwamp = gwmag * exp((r(i)-Re)/(2*H)) / exp(hmax/(2*H));
            if gwamp > gwmag, gwamp = gwmag; end
            if gwamp < -gwmag, gwamp = -gwmag; end
            
            fac = gwamp * cos(gwkh*j*dth*Re + gwkp*k*dph*Re + gwkr*(r(i)-Re));

            nd3(i,j,k) = ndt(i) * ( 1 + fac );
            
        end
    end
end


%%

h1 = figure(1);
ax1 = subplot(221);
imagesc(th*Re/1e3,(r-Re)/1e3,log10(abs(squeeze(nd3(:,:,1)))));
axis xy;

ax2 = subplot(222);
imagesc(ph*Re/1e3,(r-Re)/1e3,log10(abs(squeeze(nd3(:,1,:)))));
axis xy;

ax3 = subplot(223);
imagesc(th*Re/1e3,ph*Re/1e3,squeeze(nd3(110,:,:)));
colorbar;