function s = get3drunparams(datadir)

re = 6370e3;

% read parameters from output sferic file

fid = fopen([datadir 'sferic.dat'],'r');
s.tsteps = fread(fid,1,'int');
s.rr = fread(fid,1,'int');
s.hh = fread(fid,1,'int');
s.pp = fread(fid,1,'int');
s.numfiles = fread(fid,1,'int');
s.dt = fread(fid,1,'double');
s.r = fread(fid,s.rr,'double');
s.th = fread(fid,s.hh,'double');
s.ph = fread(fid,s.pp,'double');
s.decfactor = fread(fid,1,'int');
fclose(fid);

s.dr = diff(s.r);
s.dth = s.th(2)-s.th(1);

% repeat for inputs.dat file, and add to s structure

fid = fopen([datadir 'inputs.dat'],'r');
s.RE = fread(fid,1,'double');
s.dopml_top = fread(fid,1,'int');
s.dopml_wall = fread(fid,1,'int');
s.doionosphere = fread(fid,1,'int');
s.doioniz = fread(fid,1,'int');
s.doelve = fread(fid,1,'int');
s.dodetach = fread(fid,1,'int');
s.dotransmitter = fread(fid,1,'int');
s.savefields = fread(fid,6,'int');
s.pecground = fread(fid,1,'int');
s.maxalt = fread(fid,1,'double');
s.stepalt = fread(fid,1,'double');
s.dr0 = fread(fid,1,'double');
s.dr1 = fread(fid,1,'double');
s.dr2 = fread(fid,1,'double');
s.nground = fread(fid,1,'int');
s.range = fread(fid,1,'double');
s.drange = fread(fid,1,'double');
s.dt = fread(fid,1,'double');
s.tsteps = fread(fid,1,'int');
s.sig = fread(fid,1,'double');
s.sigm = fread(fid,1,'double');
s.camdist = fread(fid,1,'double');
s.camalt = fread(fid,1,'double');
s.elvesteps = fread(fid,1,'int');
s.sourcedirection = fread(fid,1,'int');
s.numfiles = fread(fid,1,'int');
s.venus = fread(fid,1,'int');
s.decfactor = fread(fid,1,'int');
s.nprobes = fread(fid,1,'int');
s.prober = fread(fid,s.nprobes,'int');
s.probet = fread(fid,s.nprobes,'int');
s.probep = fread(fid,s.nprobes,'int');
fclose(fid);

if s.decfactor == 0,
    s.decfactor = 1;
end

% read magnetic field components

fid = fopen([datadir 'B0.dat'],'r');
s.Br = fread(fid,1,'double');
s.Bt = fread(fid,1,'double');
s.Bp = fread(fid,1,'double');
fclose(fid);