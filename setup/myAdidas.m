function myAdidas(tag, eplv, N, nMonte)
% runDMC(tag, eplv, N, nMonte)
strtm = tic;

% create a file name == job name
bsnm = ['wrap_t' num2str(tag) '_e' num2str(round(100*eplv)) '_N' ...
        num2str(N) '_M' num2str(nMonte) ];

% open the file
fid = fopen(bsnm, 'w');
fprintf(fid, '#!/bin/bash\n');
% change to the correct directory:
fprintf(fid, 'cd /shared/users/dstrauss/stats330/hw4/part_2/\n');
% tell matlab to go 
fprintf(fid, ['matlab -nodisplay -r "rapperb(%d, %d,%d,%d); exit" '], ...
        tag,eplv,N,nMonte);
% close the file/script
fclose(fid);
%make it executable
system(['chmod +x ' bsnm]);

% create the submit command as you wish
cmd = ['qsub ' bsnm ' -l nodes=1:ppn=2 -l walltime=24:00:00'];

% submit it
[blago jbnm] = system(cmd);
jbnm(length(jbnm)) = [];
% disp(jbnm, bsnm)

% wait for the job to finish
waitforexit(jbnm, bsnm)

% control returned after the job has completed

disp(['T' num2str(tag) ' e' num2str(eplv) ' N' num2str(N) ' M' ...
      num2str(nMonte) ' time = ' num2str(toc(strtm))])


