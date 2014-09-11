function s = get1drunparams(datadir,datatype)

% read parameters from output sferic file

fid = fopen([datadir '/sferic.dat'],'r');
s.tsteps = fread(fid,1,'int');
s.xx = fread(fid,1,'int');
s.numfiles = fread(fid,1,'int');
s.dt = fread(fid,1,datatype);
s.x = fread(fid,s.xx,datatype);
fclose(fid);

s.dx = diff(s.x);

% repeat for inputs.dat file, and add to s structure

fid = fopen([datadir '/inputs.dat'],'r');
s.doionosphere = fread(fid,1,'int');
s.doioniz = fread(fid,1,'int');
s.maxalt = fread(fid,1,datatype);
s.stepalt = fread(fid,1,datatype);
s.dx1 = fread(fid,1,datatype);
s.dx2 = fread(fid,1,datatype);
s.dt = fread(fid,1,datatype);
s.tsteps = fread(fid,1,'int');
s.Esource = fread(fid,s.tsteps,datatype);
s.numfiles = fread(fid,1,'int');
fclose(fid);

fid = fopen([datadir '/B0.dat'],'r');
s.Br = fread(fid,1,datatype);
s.Bt = fread(fid,1,datatype);
s.Bp = fread(fid,1,datatype);
fclose(fid);

fid = fopen([datadir '/ne.dat'],'r');
s.ne = fread(fid,s.xx,'double');
fclose(fid);

fid = fopen([datadir '/nd.dat'],'r');
s.nd = fread(fid,s.xx,'double');
fclose(fid);

fid = fopen([datadir '/ni.dat'],'r');
s.ni = fread(fid,s.xx,'double');
fclose(fid);
