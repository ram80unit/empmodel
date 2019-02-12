% write pbs file for a given runname

function pbsfile = writepbsfile(rundir,runname,exefile)

pbsfile = [rundir '/' runname '.pbs'];

fid = fopen(pbsfile,'wt');

fprintf(fid,'#!/bin/sh\n\n');
fprintf(fid,'rm -f %s/output_K.dat\n',rundir);
fprintf(fid,'rm -f %s/output_E.dat\n',rundir);
fprintf(fid,'rm -f %s/output_D.dat\n',rundir);
fprintf(fid,'rm -f %s/output_H.dat\n',rundir);
fprintf(fid,'rm -f %s/output_J.dat\n',rundir);
fprintf(fid,'rm -f %s/output_T.dat\n',rundir);
fprintf(fid,'rm -f %s/output_O.dat\n',rundir);
fprintf(fid,'rm -f %s/output_S.dat\n',rundir);
fprintf(fid,'rm -f %s/Probe.dat\n',rundir);
fprintf(fid,'rm -f %s/elve.dat\n',rundir);
fprintf(fid,'rm -f %s/sferic.dat\n',rundir);
fprintf(fid,'\n');
fprintf(fid,'time %s/%s\n',rundir,exefile);
fclose(fid);

