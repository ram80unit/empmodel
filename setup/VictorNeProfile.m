function ne = VictorNeProfile(alt,type)

% type should be 1, 2, or 3.

filename = ['vic.amb.' num2str(type)];

A = textread(filename);

tempalts = A(1:end-1,1);
tempdens = A(1:end-1,2) * 1e6;

ne = interp1(tempalts,tempdens,alt,'cubic');

if max(alt) > max(tempalts),
    
    ind = find(alt > max(tempalts),1,'first');
    ne2 = IRIionosphere1(alt(ind:end));
    
    % multiply ne2 by constant to match ne at base altitude
    ne2 = ne2 * (ne(ind) / ne2(1));
    
    ne = [ne(1:ind-1); ne2];
end
