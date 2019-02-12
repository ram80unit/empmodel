%% launch batch of EMP 2D jobs. Use also for single jobs.

clear all; close all;

loadconstants;

% set defaults
emp2Ddefaults_disp;

% anything you want to change?
inputs.submitjob = 1;
inputs.I0 = 10e3;

dr1 = [250 500 1000];
dr2 = [250 250 500];

% master directory for set of runs
toprundir = '/projects/roma8490/emp/runs/disptest1000/';

numsims = 10;

% randomize ionosphere and magnetic field

pd = makedist('uniform','lower',72,'upper',87);
hprime = random(pd,numsims,1);
pd = makedist('uniform','lower',0.3,'upper',0.7);
beta = random(pd,numsims,1);

pd = makedist('uniform','lower',40000,'upper',50000);
B0 = random(pd,numsims,1);

% submit jobs

for m = 1:numsims,
    
    inputs.beta = beta(m);
    inputs.hprime = hprime(m);
    inputs.B0 = B0(m) * 1e-9;
    
    for k = 1:3,
	  % repeat each sim three times, for different grid sizes

	  inputs.dr1 = dr1(k);
          inputs.dr2 = dr2(k);
          inputs.drange = inputs.dr1;

          inputs.runname = sprintf('random_%03d_%d',m,k);
    
          inputs.rundir = [toprundir inputs.runname];
    
          % launch job
          [in,jobid] = emp2Drun(inputs);
    
          drawnow;

    end    
end

