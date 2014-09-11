function ne = YukiIonosphere(alt,beta,hk)

% beta and hk are Wait and Spies profile parameters
% will use IRI profile above crossover point

%alt = 0:1:120;
%hk = 90;
%beta = 0.8;
thresh = 2e9;

ne = 1.43e13 * exp(-0.15*hk) * exp((beta-0.15)*(alt-hk));

ind = find(ne > thresh,1,'first');


ne(ind:end) = thresh;
ne = smooth(ne,5);
    
%figure;
%semilogx(ne,alt,'b-');
%hold on;

%set(gca,'ylim',[70 100],'xlim',[1e6 1e10]);