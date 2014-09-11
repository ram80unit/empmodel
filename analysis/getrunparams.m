function s = getrunparams(datadir,datatype)

re = 6370e3;

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
s.dopml_top = fread(fid,1,'int');
s.dopml_wall = fread(fid,1,'int');
s.doioniz = fread(fid,1,'int');
s.doelve = fread(fid,1,'int');
s.dodetach = fread(fid,1,'int');
s.dotransmitter = fread(fid,1,'int');
s.pecground = fread(fid,1,'int');
s.maxalt = fread(fid,1,datatype);
s.stepalt = fread(fid,1,datatype);
s.dr0 = fread(fid,1,datatype);
s.dr1 = fread(fid,1,datatype);
s.dr2 = fread(fid,1,datatype);
s.nground = fread(fid,1,'int');
s.range = fread(fid,1,datatype);
s.dt = fread(fid,1,datatype);
s.tsteps = fread(fid,1,'int');
s.sig = fread(fid,1,datatype);
s.sigm = fread(fid,1,datatype);
s.camdist = fread(fid,1,datatype);
s.camalt = fread(fid,1,datatype);
s.Isource = fread(fid,s.tsteps,datatype);
s.txf0 = fread(fid,1,datatype);
s.sourcealt = fread(fid,1,datatype);
s.rsspeed = fread(fid,1,datatype);
s.numfiles = fread(fid,1,'int');
s.venus = fread(fid,1,'int');
s.decfactor = fread(fid,1,'int');
fclose(fid);