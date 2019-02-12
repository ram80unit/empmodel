function [tout,zout,Icout] = interpolate_current(t,z,Ic,zi,ti)
%
% [tout,zout,Icout] = interpolate_current(t,z,Ic,zi,ti)
%
% Function to interpolate the source current 
% according the required resolution
% defined by the vectors zi and ti
%
% tout  = ti;
% zout  = zi;
% [t1,z1] = meshgrid(t,z);
% [ti1,zi1] = meshgrid(ti,zi);
% Icout = interp2(t1,z1,Ic,ti1,zi1,'cubic');
%
% Author: Caitano L. da Silva (Penn State University)
%
% Created: Jan/30/2015
% Last modification: Jan/30/2015
%

tout  = ti;
zout  = zi;
[t1,z1] = meshgrid(t,z);
[ti1,zi1] = meshgrid(ti,zi);
Icout = interp2(t1,z1,Ic,ti1,zi1,'cubic');

return


