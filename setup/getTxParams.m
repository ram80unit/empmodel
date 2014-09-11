function in = getTxParams(in,txname)

% get all the important transmitter parameters and store them as the useful
% variables for the EMP code.

loadconstants;

txinfo = GetTxInfo(txname);

% location

in.Trlat = txinfo.lat;
in.Trlon = txinfo.lon;

% power, frequency

Prad = txinfo.power;
in.txf0 = txinfo.freq;

% convert to current

%correction April 25, 2014: 24 instead of 48 for monopole.
% another correction May 3: factor of 2 in sourcealt, since it is the
% monopole length
in.I0 = sqrt(24*pi*e0*vp^3*Prad) / (2*pi*in.txf0*(2*in.sourcealt));

% now get magnetic field direction and magnitude from IGRF

[Bx,By,Bz] = igrf(today,in.Trlat,in.Trlon,100,'geodetic');

Bh = sqrt(Bx.^2 + By.^2);

% input Bmag should be Br (vertical,upward), Bt (zero), and Bp (horizontal,
% positive)
% convert from nT to T!!
in.Bmag = [-Bz 0 Bh] * 1e-9;
