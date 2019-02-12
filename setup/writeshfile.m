% write pbs file for a given runname

function shfile = writeshfile(in)

shfile = [in.rundir '/' in.runname '.sh'];

fid = fopen(shfile,'wt');

fprintf(fid,'#!/bin/sh\n\n');
fprintf(fid,'\n');
fprintf(fid,'#SBATCH --job-name %s\n',in.runname);
fprintf(fid,'#SBATCH --qos %s\n',in.cluster);
fprintf(fid,'#SBATCH --nodes 1\n');
fprintf(fid,'#SBATCH --ntasks-per-node %s\n',in.numnodes);
fprintf(fid,'#SBATCH --time %s\n',in.walltime);
fprintf(fid,'\n');
fprintf(fid,'module load slurm\n');
fprintf(fid,'\n');
fprintf(fid,'srun %s/%s\n',in.rundir,in.exefile);
fclose(fid);

