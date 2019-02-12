% EMP ground parameters, for 2D transmitter propagation version

function [sigma,epsilon] = getGroundParams(in)

doplot = 0;

% startLat = 40;
% startLon = -100;
% az = 90;       % south
% range = 20000e3;
% dr1 = 1e3;

RE = 6370e3;  % need this for proper ground parameters along real-Earth path

npts = round(in.range/in.dr1) + 1;
interp = 0;

[stopLat,stopLon] = reckon(in.Trlat,in.Trlon,in.range*180/pi/RE,in.az);

load ( 'conductivity_data.mat' );

% compute the great circle path, assume spherical for now
[gc_lat,gc_lon] = track2(in.Trlat,in.Trlon,stopLat,stopLon,[RE 0],'degrees',npts);

% Find the distance.
dist = pi/180 * RE * distance(in.Trlat,in.Trlon,stopLat,stopLon);

% Extract the conductivities along the path
sigma = zeros(length(gc_lat),1);
epsilon = zeros(length(gc_lat),1);
nlon = length(lon);
rlon = lon(end)-lon(1);
nlat = length(lat);
rlat = lat(end)-lat(1);

for i=1:length(gc_lat)
    sigma(i) = sigmamap(floor((nlat/rlat)*(gc_lat(i)-lat(1)))+1,...
        floor((nlon/rlon)*(gc_lon(i)-lon(1)))+1);
    epsilon(i) = epsilonmap(floor((nlat/rlat)*(gc_lat(i)-lat(1)))+1,...
        floor((nlon/rlon)*(gc_lon(i)-lon(1)))+1);
end

if doplot,
    
    h1 = figure(1);
    set(h1,'position',[300 150 1000 600])
    
    ax1 = subplot(221);
    imagesc([-180:0.5:180],[-90:0.5:90],log10(sigmamap));
    axis xy;
    hold on;
    plot(gc_lon,gc_lat,'w');
    
    ax2 = subplot(223);
    plot((1:npts)*in.dr1/1e3,log10(sigma));
    
    ax3 = subplot(222);
    imagesc([-180:0.5:180],[-90:0.5:90],epsilonmap);
    axis xy;
    hold on;
    plot(gc_lon,gc_lat,'w');
    
    ax2 = subplot(224);
    plot((1:npts)*in.dr1/1e3,epsilon);
    
end

