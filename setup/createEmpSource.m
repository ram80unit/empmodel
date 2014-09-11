
function Jsv = createEmpSource(in)

testing = 0;

%in.lightningtype = 1;

decaytype = in.decaytype;

% initialize time-dependent source

Iin = zeros(in.tsteps,1);

% okay, set up current waveform: linear rise, exp decay
for t = 1:in.tsteps,
    if (t * in.dt < in.taur),
        Iin(t) = in.I0 * t * in.dt / in.taur;
    else
        Iin(t) = in.I0 * exp(-(t*in.dt-in.taur)^2/in.tauf^2) + in.Ic*(1 - exp(-(t*in.dt-in.taur)^2/in.tauf^2));
    end
end

% okay, filter it

b = fir1(40,in.fcut/(1/in.dt/2));

Iin2 = filter(b,1,Iin);
Iin2 = Iin2 / max(Iin2) * in.I0;   % rescaled to peak current

% fix delays and stuff
if decaytype < 3,
    in.Iin = [zeros(round(in.sourcealt/abs(in.rsspeed)/in.dt),1); Iin2];
else
    in.Iin = [Iin2; zeros(2*round(in.sourcealt/abs(in.rsspeed)/in.dt),1)];
end

in.Iin(in.Iin < 1e-6) = 1e-6;

in.tsteps2 = length(in.Iin);


%% if transmitter

if in.dotransmitter,
    decaytype = 7;
    in.sourcealt = 4e3;
    in.lightningtype = 0;
    in.Iin = in.I0 * sin(2*pi*in.txf0 * in.dt * (1:in.tsteps));
    in.tsteps2 = in.tsteps;
    
end


%% spatial variation

% find index of top of source
if in.lightningtype == 0,
    satop = in.nground + floor(in.sourcealt/in.dr1) + 1;
else
    satop = in.nground + floor(in.sourcealt/in.dr1 + in.chlength/2/in.dr1) + 1;
end

% a few cells above and below the top
nextra = 0;

% total length to top of defined source. Will write this to source file
nstotal = satop + nextra;

% setup output array
Jsv = zeros(nstotal,in.tsteps2);

% create r vector
r = zeros(nstotal,1);
r(1) = in.Re - in.nground*in.dr0;
for m = 2:in.nground+1,
    r(m) = r(m-1) + in.dr0;
end
for m = in.nground+2:nstotal,
    r(m) = r(m-1) + in.dr1;
end


%% CG source
if in.lightningtype == 0,
    
    % altitude-dependent delay of return stroke from initiation
    rsdelay = round((r-in.Re)/abs(in.rsspeed)/in.dt);
    
    for t = 1:in.tsteps2,
        
        for i = 1:nstotal,
            % index into source delay
            tJs = t + max(rsdelay) - rsdelay(i);
            if tJs > in.tsteps2,
                tJs = in.tsteps2;
            end
            
            if decaytype == 0,       % TL source
                
                Jsv(i,t) = in.Iin(tJs);
                
            elseif decaytype == 1,   % MTLL
                
                Jsv(i,t) = in.Iin(tJs) * (1 - (r(i)-in.Re)/in.sourcealt);
                
            elseif decaytype == 2,   % MTLE
                
                rssigma = in.sourcealt/3;      % decay length for MTLE
                Jsv(i,t) = in.Iin(tJs) * exp(-(r(i)-in.Re)/rssigma);
                
            elseif decaytype == 3,   % BG
                
                if t > rsdelay(i),
                    Jsv(i,t) = in.Iin(t);
                end
                
            elseif decaytype == 4,   % TCS
                
                tcsdelay = round((r(i)-in.Re)/3e8/in.dt);
                tJs = t + tcsdelay;
                if tJs > in.tsteps2, tJs = in.tsteps2; end
                if t > rsdelay(i),
                    Jsv(i,t) = in.Iin(tJs);
                end
                
            elseif decaytype == 5,   % DU
                
                taud = 2e-6;        % decay constant for DU source
                tcsdelay = round((r(i)-in.Re)/3e8/in.dt);
                tJs = t + tcsdelay;
                vstar = in.rsspeed/(1+abs(in.rsspeed)/3e8);
                z = (i-1)*in.dr1;
                tcsind = round(z/vstar/in.dt);
                if tcsind < 1, tcsind = 1; end
                if tJs > in.tsteps2, tJs = in.tsteps2; end
                if t > rsdelay(i),
                    Jsv(i,t) = (in.Iin(tJs) - exp(-(t*in.dt-z/abs(in.rsspeed))/taud)*in.Iin(tcsind));
                end
                
            elseif decaytype == 6,   % dummy
                
                Jsv(i,t) = in.Iin(t);
                
            elseif decaytype == 7,   % transmitter: decay with alt, no propagation
                
                Jsv(i,t) = in.Iin(t) * (1 - (r(i)-in.Re)/in.sourcealt);
                
            end
            
        end
        
    end
    
    %% IC source
elseif in.lightningtype == 1,
    
    % for IC, we will assume current travels from one end to the other with
    % no reflection, and no decay - simply disappears at far end.
    
    basealt = in.sourcealt - in.chlength/2;
    topalt = in.sourcealt + in.chlength/2;
    
    i1 = find((r-in.Re) >= basealt, 1, 'first');
    i2 = find((r-in.Re) >= topalt, 1, 'first');
    
    % altitude-dependent delay of return stroke from initiation
    rsdelay = round((r-in.Re)/abs(in.rsspeed)/in.dt);
    
    for t = 1:in.tsteps2,
        
        for i = i1:i2,
            % index into source delay
            tJs = t + max(rsdelay) - rsdelay(i-i1+4);
            if tJs > in.tsteps2,
                tJs = in.tsteps2;
            end
            Jsv(i,t) = in.Iin(tJs) * (i2-i)/(i2-i1);        % MTLL-like decay downwards            
        end
        
    end
    
    % if rsspeed is negative, downwards, so flip the good part
    if in.rsspeed < 0,
        Jsv(i1:i2,:) = flipud(Jsv(i1:i2,:));
    end
    
    
    %% CID source
elseif in.lightningtype == 2,
    
    
    
end


%% we only need to keep Jsv up to about 0.1 ms. Find it automatically though.

ind = find(mean(Jsv,1) > 1e-6, 1, 'last');
Jsv = Jsv(:,1:ind);

if testing,
    figure;
    imagesc(in.dt*(1:ind)*1e3,in.dr1*(1:nstotal)/1e3,Jsv);
    axis xy;
end

