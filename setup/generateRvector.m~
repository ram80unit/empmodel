function [r,dr] = generateRvector(dr0,dr1,dr2,nground,stepalt,maxalt)

loadconstants;

rr = stepalt/dr1 + (maxalt - stepalt)/dr2 + 1 + nground;
r = zeros(rr,1);
dr = zeros(rr-1,1);
r(1) = RE - nground*dr0;
for i = 2:nground+1,
    r(i) = r(i-1) + dr0;
end
for i = nground+2:rr,
    if r(i-1) < (RE + stepalt),
        r(i) = r(i-1) + dr1;
    else
        r(i) = r(i-1) + dr2;
    end
    dr(i-1) = r(i) - r(i-1);
end