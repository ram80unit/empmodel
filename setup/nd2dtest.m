% 2D nd test

clear all; close all;

Re = 6370e3;

r = Re:500:(Re+110e3);
range = 0:500:300e3;
th = range/Re;
dth = th(2) - h(1);
H = 6e3;
hmax = 100e3;
gwmag = 0.5;
gwlam = 20e3;

nd = MSISatmosphere1((r-Re)/1000);
ndt = nd.total * 1e6;

nd2 = zeros(length(r),length(h));

for i = 1:length(r),
    for j = 1:length(h),
        
        fac = gwmag / exp(hmax/(2*H)) * exp((r(i)-Re)/(2*H)) * sin(2*pi/gwlam*j*dh*Re);
        if fac > gwmag, fac = gwmag; end
        if fac < -gwmag, fac = -gwmag; end
        nd2(i,j) = ndt(i) * ( 1 + fac );

    end
end

imagesc(range/1000,(r-Re)/1000,log10(abs(nd2)));
axis xy;
colorbar;
colormap(jet(128));