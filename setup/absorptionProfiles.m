% calculate and plot absorption curves (nu * ne)

clear all; close all;

% do a bunch of different ionospheres, but just one atmosphere

QE = 1.6e-19;
ME = 9.11e-31;

alt = 0:1:150;

beta = [0.5 0.7 0.9 1.1];
hk = [84:1:88];

ne0 = IRIionosphere1(alt);

relfac = [0.1 0.2 0.5 1 2];

ne = zeros(length(alt),length(relfac)+length(beta)*length(hk));

for m = 1:length(relfac),
    ne(:,m) = ne0 * relfac(m);
end


for i = 1:length(beta),
    for j = 1:length(hk),
        m = length(relfac) + j + (i-1)*length(hk);
        ne(:,m) = YukiIonosphere(alt,beta(i),hk(j));
    end
end

nd = MSISatmosphere1(alt);
ndt = nd.total * 1e6;
mue = 1.4856 * ndt(1) ./ ndt;
nu = QE/ME ./ mue;


absorption = zeros(size(ne));

for m = 1:size(ne,2),
    absorption(:,m) = nu' .* ne(:,m);
end

figure;
plot(absorption,alt);
xlabel('Absorption (a.u.)');
ylabel('Altitude (km)');

ylim([50 120]);
