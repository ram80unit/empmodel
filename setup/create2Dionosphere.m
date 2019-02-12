
function in = create2Dionosphere(in)

loadconstants;

% read UF h', beta values and create ne array in 2D.

% need to create distance array
thmax = in.range / in.Re;
dth = in.drange / in.Re;
hh = round(thmax / dth) + 1;
th = 0:dth:(hh-1)*dth;
dvec = th*in.Re;

% test for comparison
dvec2 = 0:in.drange:in.range;

%% Rasoul perturbations, given by files

if strcmp(in.iono2Dmethod,'eclipse'),

    A = load(in.iono2Dfile);
    in.ne = squeeze(A.CO.e(in.eclipseT, :, :));
    
    mue = 1.4856 * in.ndt(in.nground+1) ./ in.ndt;
    nue = (QE / ME) ./ mue;
    for i = 1:length(in.ne(:,1))
        in.nu(i, 1:length(in.ne(1, :))) = nue;
    end

    in.ne = in.ne';
    in.nu = in.nu';
    
%% Precipitating perturbations, given by 1D profile, but also given a width as input.

elseif strcmp(in.iono2Dmethod,'Precip'),

    % read and interpolate 2D ionosphere from Rasoul runs
    
    A = load(in.iono2Dfile);
    % gives alt (km), Ne0 (per cc), eprofs (nalt x 5, per cc)
    
    % use Ne0 to create basic 2D ionosphere
    rout = (in.r-in.Re)/1e3;
    ne0 = interp1(A.alt0,A.Ne0'*1e6,rout,'pchip');
    
    ne02d = repmat(ne0,1,in.hh);
    
    % now add perturbation
    
    nepert0 = (A.eprof(:,in.iono2Dindex) - A.Ne0')*1e6;
    nepert = interp1(A.alt0,nepert0,rout,'pchip');
    
    % need distance of each of hh grid cells along ground
    dist = in.drange * (0:1:(in.hh-1));
    
    for m = 1:in.hh, 
        thispert = nepert * exp(-(dist(m) - in.ionopertdist)^2 / in.ionopertsize^2);
        ne02d(:,m) = ne02d(:,m) + thispert;
    end
    
    
    in.ne = ne02d;
    
    
    
    
    %% similar for beams, but 1D profiles in one grid cell.
    
elseif strcmp(in.iono2Dmethod,'beam'),
    
    A = load(in.iono2Dfile);
    % gives alt,nebg,nepert.
    % cut it off above 80 km, to maximum of 3e3.
    for k = 81:101,
        if (A.nebg(k) > 3e3), 
            A.nebg(k) = 3e3;
        end
        if (A.nepert(k) > 3e3),
            A.nepert(k) = 3e3;
        end
    end
    
    rout = (in.r-in.Re)/1e3;
    ne0 = interp1(A.alt,A.nebg*1e6,rout);

    ne0(isnan(ne0)) = max(A.nebg(:))*1e6;
    ne02d = repmat(ne0,1,hh);

%if ~strfind(in.runname,'amb'),    
   
 i1 = find(dvec >= in.ionopertdist, 1, 'first');
disp(i1);
    
    nepert2 = interp1(A.alt,A.nepert*1e6,rout);
    
    nepert2(isnan(nepert2)) = max(A.nebg(:))*1e6;

    ne02d(:,i1) = nepert2;
    ne02d(:,i1-1) = ne0 + (nepert2-ne0)*0.3;
    ne02d(:,i1+1) = ne0 + (nepert2-ne0)*0.3;
    
%end
    in.ne = ne02d;
    
    
    %% Basic variation in h', beta. 
    % I have inputs called in.hbrange, in.hparray, in.betaarray.
    
else
    % initialize
    ne0 = YukiIonosphere((in.r-in.Re)/1e3,in.betaarray(1),in.hparray(1));
    ne02d = repmat(ne0,1,hh);
    
    for m = 2:length(in.hbrange),
        ne1 = YukiIonosphere((in.r-in.Re)/1e3,in.betaarray(m),in.hparray(m));
        ind = find(dvec > in.hbrange(m),1,'first');
        if ~isempty(ind),
            ne02d(:,ind:end) = repmat(ne1,1,(hh-ind+1));
        end
    end
    
    in.ne = ne02d;
    
end

