%% launch batch of EMP 3D jobs. Use also for single jobs.

clear all; close all;

loadconstants;

% set defaults
emp3DdefaultsTx;

% anything you want to change?
inputs.submitjob = 1;
inputs.doioniz = 0;
inputs.dodetach = 0;
inputs.doionosphere = 0;
Prad = 1e6;
inputs.cluster = 'batchnew';
inputs.numnodes = '12';

% master directory for set of runs
toprundir = '/shared/users/ram80/empcodes/runs/transtest/';

% variable for batch of runs. name must match an input!
var1.name = 'TX';
%var1.values = [{'NPM'} {'NAA'} {'NAU'} {'NWC'} {'NLK'}];
var1.values = [{'NAA'}];


%% submit jobs

for m = 1:length(var1.values),
    
    if isa(var1.values,'cell'),
       
        % get transmitter power->current, lat, lon, dip angle, Bmag,
        % frequency
        if strcmp(var1.name,'TX'),
            inputs = getTxParams(inputs,var1.values(m));
        end
        
        inputs.runname = [var1.name '_' var1.values{m}];
        
    elseif isa(var1.values(m),'double'),
        
        evalstr = ['inputs.' var1.name ' = ' num2str(var1.values(m)) ';'];
        eval(evalstr);
        
        if exist('Prad','var'),
            inputs.I0 = sqrt(24*pi*e0*vp^3*Prad) / (2*pi*inputs.txf0*(2*inputs.sourcealt));
        end  
        
        inputs.runname = [var1.name '_' sprintf('%03d',var1.values(m))];  
        
    end
    
    inputs.runname(strfind(inputs.runname,'+')) = '';
    inputs.rundir = [toprundir inputs.runname];
    
    % launch job
    [in,jobid] = emp3Drun(inputs);
    
end
