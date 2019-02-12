% compare ionosphere profiles

clear all; close all;

loadconstants;

r = 0:200:110e3;
beta = 0.5:0.2:1.1;
hk = 87;

ne1 = IRIionosphere1(r/1000);

nd = MSISatmosphere1(r/1000);
ndt = nd.total * 1e6;

f = 20e3;

mu = 1.48 * ndt(1)./ndt;
nu = QE ./ (mu*ME);

wp1sq = QE^2 * ne1 / (ME * e0);

h1 = figure(1);
ax = axes;

semilogx(ax,ne1,r/1e3,'b');
hold(ax,'on');
semilogx(ax,nu,r/1e3,'k');

A = wp1sq./(nu*2*pi*f);
i1 = find(A > 1, 1, 'first');

plot(ax,[1e0 1e10],[r(i1) r(i1)]/1e3,'b:');


%%

for m = 1:length(beta),
    
ne2 = YukiIonosphere(r/1000,beta(m),hk)';

wp2sq = QE^2 * ne2 / (ME * e0);

semilogx(ax,ne2,r/1e3,'r');

A = wp2sq./(nu*2*pi*f);
i1 = find(A > 1, 1, 'first');

plot(ax,[1e0 1e10],[r(i1) r(i1)]/1e3,'r:');

end

set(ax,'xlim',[1e0 1e10],'ylim',[40 105]);

