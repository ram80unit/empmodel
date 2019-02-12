% compare ionospheres

clear all; close all;

RE = 6370e3;

dr = 1000;
r = RE:dr:(RE+150e3);

h1 = figure(1); set(h1,'position',[200 200 600 600]);
ax = axes;

ne1 = IRIionosphere1((r-RE)/1000);
ne2 = IRIionosphere2((r-RE)/1000);
ne3 = IRIionosphere3((r-RE)/1000);

ne4 = ne1 * 0.7;


semilogx(ne1,(r-RE)/1000,'k');
hold(ax,'on');
semilogx(ne2,(r-RE)/1000,'b');
semilogx(ne3,(r-RE)/1000,'r');
semilogx(ne4,(r-RE)/1000,'k:');

xlabel(ax,'Electron density per m3');
ylabel(ax,'Altitude (km)');
legend(ax,'IRI 1','IRI 2','IRI 3','half IRI 1');