function rates = getNonlinearRates(r,ndcm,nground)

% r input is altitude in km. Not to be confused with alt below, which is
% the hard coded altitudes in bolsig+ outputs.

doplot = 0;

alt = [-25:25:400]';

nd = ndcm.total * 1e6;
N0 = nd(nground+1);        

ioniz = zeros(length(alt),100);
attach = zeros(length(alt),100);
mobility = zeros(length(alt),100);
energy = zeros(length(alt),100);
Ored = zeros(length(alt),100);
Ogrn = zeros(length(alt),100);
N21p = zeros(length(alt),100);
N22p = zeros(length(alt),100);
N2p1N = zeros(length(alt),100);
O2p1N = zeros(length(alt),100);
N2pM = zeros(length(alt),100);

for m = 1:length(alt),
    
    k = find(r <= alt(m),1,'last');
    if isempty(k),
        k = 1;
    end
    
    reldens.o = ndcm.o(k) ./ ndcm.total(k);
    reldens.n2 = ndcm.n2(k) ./ ndcm.total(k);
    reldens.o2 = ndcm.o2(k) ./ ndcm.total(k);
    reldens.ar = ndcm.ar(k) ./ ndcm.total(k);
    reldens.he = ndcm.he(k) ./ ndcm.total(k);
    
    bolsigfile = sprintf('air_%03dkm.dat',alt(m));
    
    [ea, air] = readoutputs_air(bolsigfile,doplot,reldens);
    
    ioniz(m,:) = air.ioniz;
    attach(m,:) = air.attach;
    mobility(m,:) = air.mobility;
    energy(m,:) = air.energy;
    Ored(m,:) = air.Ored;
    Ogrn(m,:) = air.Ogrn;
    N21p(m,:) = air.n21p;
    N22p(m,:) = air.n22p;
    N2p1N(m,:) = air.mossn2p1n;
    O2p1N(m,:) = air.mosso2p1n;
    N2pM(m,:) = air.mossn2pm;
    
end

% threshold for mobility
for m = 1:length(alt),
    for n = 1:length(ea),
        if mobility(m,n) > 1.4856 * N0,
            mobility(m,n) = 1.4856 * N0;
        end
    end
end

% interpolate those 2D arrays on useful r vector

ind = find(ea*N0 > 1e7,1,'first');

rates.ioniz = interp2(ea,alt,ioniz,ea(1:ind),r,'linear');
rates.attach = interp2(ea,alt,attach,ea(1:ind),r,'linear');
rates.mobility = interp2(ea,alt,mobility,ea(1:ind),r,'linear');
rates.energy = interp2(ea,alt,energy,ea(1:ind),r,'linear');
rates.Ored = interp2(ea,alt,Ored,ea(1:ind),r,'linear');
rates.Ogrn = interp2(ea,alt,Ogrn,ea(1:ind),r,'linear');
rates.N21p = interp2(ea,alt,N21p,ea(1:ind),r,'linear');
rates.N22p = interp2(ea,alt,N22p,ea(1:ind),r,'linear');
rates.N2p1N = interp2(ea,alt,N2p1N,ea(1:ind),r,'linear');
rates.O2p1N = interp2(ea,alt,O2p1N,ea(1:ind),r,'linear');
rates.N2pM = interp2(ea,alt,N2pM,ea(1:ind),r,'linear');
rates.efield = ea(1:ind);

rates.ioniz(isnan(rates.ioniz)) = 0;
rates.attach(isnan(rates.attach)) = 0;
rates.mobility(isnan(rates.mobility)) = 0;
rates.energy(isnan(rates.energy)) = 0;
rates.Ored(isnan(rates.Ored)) = 0;
rates.Ogrn(isnan(rates.Ogrn)) = 0;
rates.N21p(isnan(rates.N21p)) = 0;
rates.N22p(isnan(rates.N22p)) = 0;
rates.N2p1N(isnan(rates.N2p1N)) = 0;
rates.O2p1N(isnan(rates.O2p1N)) = 0;
rates.N2pM(isnan(rates.N2pM)) = 0;