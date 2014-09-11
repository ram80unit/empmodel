% Venus atmosphere, from Venus fact sheet

function nd = VenusAtmosphere(alt)

%Surface density:
sfden = 65.0;   % kg/m3
%Scale height: 
sch = 15.9;     % km
%Total mass of atmosphere: 
totalatm = 4.8e20;  % kg
%Average temperature: 
avtemp = 737;   % K (464 C)
%Diurnal temperature range: ~0 
%Wind speeds: 0.3 to 1.0 m/s (surface)
%Mean molecular weight: 
mmw = 43.45;    % g/mole
%Atmospheric composition (near surface, by volume): 
%    Major:       96.5% Carbon Dioxide (CO2), 3.5% Nitrogen (N2) 
%    Minor (ppm): Sulfur Dioxide (SO2) - 150; Argon (Ar) - 70; Water (H2O) - 20;
%                 Carbon Monoxide (CO) - 17; Helium (He) - 12; Neon (Ne) - 7

avogadro = 6.022e23;

% to get from surface density to surface number density, use mmw and avo

n0 = sfden * 1000 / mmw * avogadro;

n0earth = 2.688e25;

nd.total = n0 * exp(-alt/sch);
nd.temp = zeros(size(alt));

% semilogx(nd,alt);