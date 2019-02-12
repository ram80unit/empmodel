% get elevation angle to center of elve, given distance and ionosphere
% height

function elev = elevToElve(range,hiono)

Re = 6370e3;

% Earth interior angle
theta = range/Re;

d2 = sqrt(Re^2 + (Re+hiono)^2 - 2*Re*(Re+hiono)*cos(theta));

% need to get farther angle first, since it must be less than 90 deg

psi = asin(Re * sin(theta) / d2);

phi = pi - theta - psi;

% output elevation in degrees

elev = phi*180/pi - 90;