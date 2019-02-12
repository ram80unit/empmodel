function [r,dr] = generateRvector(in)

loadconstants;

rr = in.stepalt/in.dr1 + (in.maxalt - in.stepalt)/in.dr2 + 1 + in.nground;
r = zeros(rr,1);
dr = zeros(rr-1,1);
r(1) = in.Re - in.nground*in.dr0;
for i = 2:in.nground+1,
    r(i) = r(i-1) + in.dr0;
end
for i = in.nground+2:rr,
    if r(i-1) < (in.Re + in.stepalt),
        r(i) = r(i-1) + in.dr1;
    else
        r(i) = r(i-1) + in.dr2;
    end
    dr(i-1) = r(i) - r(i-1);
end
