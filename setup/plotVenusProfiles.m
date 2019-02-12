% plot Venus profiles

clear all; close all;

alt = 0:10:500;

ne = zeros(length(alt),3);

for ind = 1:3,        

ne(:,ind) = VenusIonosphere1(alt,ind);

end

nd = VenusAtmosphere(alt);

h1 = figure(1);
set(h1,'position',[1800 100 1000 600]);
ax(1) = subplot(121);
ax(2) = subplot(122);

semilogx(ax(1),ne,alt);
semilogx(ax(2),nd,alt);