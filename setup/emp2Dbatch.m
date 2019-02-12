%% launch batch of EP 2D jobs. Use also for single jobs.

clear all; close all;

loadconstants;

% set defaults
emp2Ddefaults;

% anything you want to change?
inputs.submitjob = 1;
inputs.I0 = 10e3;

% master directory for set of runs
toprundir = '/projects/wexu6668/emp/runs/eclipse';

% variable for batch of runs. name must match an input!
var1.name = 'e';
var1.values = 1:100;
inputs.eclipseTime = var1.values;

% submit jobs
for m = 1:length(var1.values),
    
    inputs.eclipseT = inputs.eclipseTime(m);
    % change variables as requested
    evalstr = ['inputs.' var1.name ' = ' num2str(var1.values(m)) ';'];
    eval(evalstr);
    
    inputs.runname = [var1.name];
    inputs.runname(strfind(inputs.runname,'+')) = '';
    sourceName = sprintf('%3.0f', var1.values(m));
    inputs.rundir = [toprundir sourceName];
    
    % launch job
    [in,jobid] = emp2Drun(inputs);
    
    drawnow;
    
end
