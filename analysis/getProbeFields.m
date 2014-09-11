% get probes

function probe = getProbeFields(datadir,s,datatype)

fid = fopen([datadir 'Probe.dat'],'r');
nprobes = fread(fid,1,'int');
probe.r = fread(fid,nprobes,'int');
probe.t = fread(fid,nprobes,'int');
% only difference between 2D and 3D is probep
if isfield(s,'probep'),
    probe.p = fread(fid,nprobes,'int');
end
probe.Er = fread(fid,[nprobes s.tsteps],datatype);
probe.Et = fread(fid,[nprobes s.tsteps],datatype);
probe.Ep = fread(fid,[nprobes s.tsteps],datatype);
probe.Hr = fread(fid,[nprobes s.tsteps],datatype);
probe.Ht = fread(fid,[nprobes s.tsteps],datatype);
probe.Hp = fread(fid,[nprobes s.tsteps],datatype);
fclose(fid);

% do magnitudes here, just for fun.
% note: some error is introduced by fields not being co-located. This would
% have to be fixed in the code

probe.Emag = sqrt(probe.Er.^2 + probe.Et.^2 + probe.Ep.^2);
probe.Hmag = sqrt(probe.Hr.^2 + probe.Ht.^2 + probe.Hp.^2);