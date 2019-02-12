% test Venus Daytime ionosphere

clear all; close all;

alt = 0:1:500;

h1 = figure(1);
set(h1,'position',[1800 100 1200 600]);

ne(:,1) = VenusIonosphere1(alt,1);

ne(:,2) = VenusIonosphere1(alt,2);

ne(:,3) = VenusIonosphere1(alt,3);

semilogx(ne,alt);