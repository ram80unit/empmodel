function s = get2drunparams(datadir,datatype)

% read parameters from output sferic file

fid = fopen([datadir '/sferic.dat'],'r');
s.tsteps = fread(fid,1,'int');
s.rr = fread(fid,1,'int');
s.hh = fread(fid,1,'int');
s.numfiles = fread(fid,1,'int');
s.dt = fread(fid,1,datatype);
s.r = fread(fid,s.rr,datatype);
s.th = fread(fid,s.hh,datatype);
s.decfactor = fread(fid,1,'int');
fclose(fid);

s.dr = diff(s.r);
s.dth = s.th(2)-s.th(1);

% repeat for inputs.dat file, and add to s structure

fid = fopen([datadir '/inputs.dat'],'r');
s.RE = fread(fid,1,datatype);
s.dopml_top = fread(fid,1,'int');
s.dopml_wall = fread(fid,1,'int');
s.doionosphere = fread(fid,1,'int');
s.doioniz = fread(fid,1,'int');
s.doelve = fread(fid,1,'int');
s.dodetach = fread(fid,1,'int');
s.dotransmitter = fread(fid,1,'int');
s.savefields = fread(fid,6,'int');
s.groundmethod = fread(fid,1,'int');
s.maxalt = fread(fid,1,datatype);
s.stepalt = fread(fid,1,datatype);
s.dr0 = fread(fid,1,datatype);
s.dr1 = fread(fid,1,datatype);
s.dr2 = fread(fid,1,datatype);
s.nground = fread(fid,1,'int');
s.range = fread(fid,1,datatype);
s.drange = fread(fid,1,datatype);
s.dt = fread(fid,1,datatype);
s.tsteps = fread(fid,1,'int');
s.sig = fread(fid,1,datatype);
s.sigm = fread(fid,1,datatype);
s.camdist = fread(fid,1,datatype);
s.camalt = fread(fid,1,datatype);
s.elvesteps = fread(fid,1,'int');
s.numfiles = fread(fid,1,'int');
%s.sourcedirection = fread(fid,1,'int');
s.planet = fread(fid,1,'int');
s.decfactor = fread(fid,1,'int');
s.nprobes = fread(fid,1,'int');
s.prober = fread(fid,s.nprobes,'int');
s.probet = fread(fid,s.nprobes,'int');
fclose(fid);

fid = fopen([datadir '/B0.dat'],'r');
s.Br = fread(fid,s.hh,datatype);
s.Bt = fread(fid,s.hh,datatype);
s.Bp = fread(fid,s.hh,datatype);
fclose(fid);

fid = fopen([datadir '/ne.dat'],'r');
s.ne = fread(fid,s.rr,'double');
fclose(fid);

fid = fopen([datadir '/nd.dat'],'r');
s.nd = fread(fid,s.rr,'double');
fclose(fid);

fid = fopen([datadir '/ni.dat'],'r');
s.ni = fread(fid,s.rr,'double');
fclose(fid);

fid = fopen([datadir '/etemp.dat'],'r');
s.Te = fread(fid,s.rr,'double');
fclose(fid);

fid = fopen([datadir '/ground.dat'],'r');
s.sigma = fread(fid,s.hh,'double');
s.epsilon = fread(fid,s.hh,'double');
fclose(fid);

if exist([datadir '/source.dat'],'file'),
    fid = fopen([datadir '/source.dat'],'r');
    s.nsalts = fread(fid,1,'int');
    s.nstimes = fread(fid,1,'int');
    s.channelcells = fread(fid,1,'int');
    s.Isource = fread(fid,[s.nsalts s.nstimes],'double');
    fclose(fid);
end

