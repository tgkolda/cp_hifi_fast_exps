%%RUN_ALIGNED_VORTEX Runs aligned vortex experiment.
% By setting the flag HASDIRECT=0, the direct methods (slow) are not evaluated
%
% See also CP_ALS_HIFI.

% Johannes Brust          johannesbrust@yahoo.com
% School of Mathematical and Statistical Sciences
% Arizona State University
%
% 01/09/26, J.B., store run information by solver

%% Check paths
% Ensure tensor_toolbox is in path
if ~exist('tensor','file')
	warning('Tensor Toolbox not found in path. Attempting to add it now...');
    addpath(genpath('../../../tensor_toolbox'));
end

% Add CP-ALS-HIFI to path
if ~exist('cp_als_hifi','file')
	warning('CP-ALS-HIFI not found in path. Adding it now...');
	addpath(genpath('../../cp_hifi'));
end

%% Setup for vortex data

% Loading data
fprintf('Loading vortex data...\n');
% Add CP-ALS-HIFI to path
if ~exist('vortex.mat','file')
	warning('Vortex data not found. Trying to add it now...');
	dpath = genpath('../../../tensor_data_vortex_shedding');
    addpath(dpath);
end

load('vortex.mat')
Xorig = tensor(X);
clear X

pname   = 'vortex';
ptype   = 'aligned';
datadir = '../results/';

%% General Setup
hifiinfo{1}.inf      = true;
hifiinfo{1}.kfunc   = kernfunc_gaussian(4); 
hifiinfo{1}.lambda  = 0.1;
hifiinfo{2}.inf      = true;
hifiinfo{2}.kfunc   = kernfunc_gaussian(4); 
hifiinfo{2}.lambda  = 0.1;  
hifiinfo{3}.inf      = true;
hifiinfo{3}.kfunc   = kernfunc_gaussian(3); 
hifiinfo{3}.lambda  = 0.1;

%% Experiment setup
rng(0,'twister');
X       = tensor_aligned(Xorig);

R       = [5,10,15,20,25,30,35,40,45,50];
nruns   = 3;
genargs = {'printitn',10,'tol',1e-6,'maxiters',50,'trace',true};
seeds   = maxk(primes(100000),nruns);

hasdirect = 0;

% Print experiment infos and start diary
rr  = matlabRelease;
ss  = strsplit(char(datetime('now')));
dn  = [ss{1},'-',strrep(ss{2},':','_')];
diary([datadir,pname,'_',ptype,'_','log','_',dn,'.txt']);

fprintf('\n*************** %s experiment **********          \n',[pname,'_',ptype]);
fprintf('*           Software: Matlab %s                     \n',rr.Release);
fprintf('*             System: %s                            \n',computer);
fprintf('*               File: %s                            \n',mfilename);
fprintf('*                                                   \n');
fprintf('*         nruns(ALS): %s                            \n',num2str(nruns));
fprintf('*              seeds: %s                            \n',num2str(seeds));
fprintf('*              ranks: %s                            \n',num2str(R));
fprintf('*                                                   \n');
fprintf('*           datetime: %s                            \n',string(datetime('now')));
fprintf('*         J.J. Brust, T.G., Kolda                   \n');
fprintf('*************************************************** \n');
fprintf('\n');

%% Run and store fast solvers
solveargsFAST    = cell(2,1);
solveargsFAST{1} = {'solver','pcg','inmaxit', 75, 'intol', 1e-6};
solveargsFAST{2} = {'solver','direct_decoupled'};

for ii=1:length(solveargsFAST)

    solveargs       = cell(1,1);
    solveargs{1}    = solveargsFAST{ii};
    
    best_exp_tab    = table('Size',[length(R) 6], ...
        'VariableTypes',{'double','double','double','double','double','double'}, ...
        'VariableNames',{'rank'  , 'runid', 'seed'  ,'time'  ,'rerr'  ,'it'});
    all_exp_tab     = table('Size',[nruns*length(R) 6], ...
        'VariableTypes',{'double','double','double','double','double','double'}, ...
        'VariableNames',{'rank'  , 'runid', 'seed'  ,'time'  ,'rerr'  ,'it'});
    
    for i=1:length(R)
        fnmstem = [datadir,pname,'_',ptype,'_'];
        [best_runs_tab, all_runs_tab] = run_aligned(X,R(i),hifiinfo,genargs,solveargs,nruns,seeds);
    
        best_exp_tab{i,1}       = R(i);
        best_exp_tab{i,2:end}   = [best_runs_tab.RunId(:) best_runs_tab.Seed(:) best_runs_tab.Time(:) best_runs_tab.RelErr(:) best_runs_tab.Iters(:)];
        
        idx                     = ((i-1)*nruns+1):(i*nruns);
        all_exp_tab{idx,1}      = R(i);
        all_exp_tab{idx,2:end}  = [all_runs_tab.RunId(:) all_runs_tab.Seed(:) all_runs_tab.Time(:) all_runs_tab.RelErr(:) all_runs_tab.Iters(:)];
    end
    
    % save tables
    fname = fullfile([datadir,pname,'_',ptype,'_',solveargs{1}{2},'_best.txt']);
    fid=fopen(fname,'w+');
    fprintf(fid,'%s \n',sprintf('%% file: %s',mfilename));
    fprintf(fid,'%s \n',sprintf('%% date: %s',string(datetime('now'))));
    writetable(best_exp_tab, fname,'WriteMode','append',"WriteVariableNames",true);
    
    fname = fullfile([datadir,pname,'_',ptype,'_',solveargs{1}{2},'_all','.txt']);
    fid=fopen(fname,'w+');
    fprintf(fid,'%s \n',sprintf('%% file: %s',mfilename));
    fprintf(fid,'%s \n',sprintf('%% date: %s',string(datetime('now'))));
    fclose(fid);
    writetable(all_exp_tab, fname, 'WriteMode','append',"WriteVariableNames",true);

end


%% Run direct solver (slow)

if hasdirect
    solveargs       = cell(1,1);
    solveargs{1}    = {'solver','direct'};

    best_exp_tab    = table('Size',[length(R) 6], ...
        'VariableTypes',{'double','double','double','double','double','double'}, ...
        'VariableNames',{'rank'  , 'runid', 'seed'  ,'time'  ,'rerr'  ,'it'});
    all_exp_tab     = table('Size',[nruns*length(R) 6], ...
        'VariableTypes',{'double','double','double','double','double','double'}, ...
        'VariableNames',{'rank'  , 'runid', 'seed'  ,'time'  ,'rerr'  ,'it'});

    runtrace         = cell(length(R),2);      

    for i=1:length(R)
        fnmstem = [datadir,pname,'_',ptype,'_'];
        [best_runs_tab, all_runs_tab, traces] = run_aligned(X,R(i),hifiinfo,genargs,solveargs,nruns,seeds);

        best_exp_tab{i,1}       = R(i);
        best_exp_tab{i,2:end}   = [best_runs_tab.RunId(:) best_runs_tab.Seed(:) best_runs_tab.Time(:) best_runs_tab.RelErr(:) best_runs_tab.Iters(:)];

        idx                     = ((i-1)*nruns+1):(i*nruns);
        all_exp_tab{idx,1}      = R(i);
        all_exp_tab{idx,2:end}  = [all_runs_tab.RunId(:) all_runs_tab.Seed(:) all_runs_tab.Time(:) all_runs_tab.RelErr(:) all_runs_tab.Iters(:)];

        [runtrace{i,:}]         = deal(R(i), traces);
    end

    % save tables
    fname = fullfile([datadir,pname,'_',ptype,'_',solveargs{1}{2},'_best.txt']);
    fid=fopen(fname,'w+');
    fprintf(fid,'%s \n',sprintf('%% file: %s',mfilename));
    fprintf(fid,'%s \n',sprintf('%% date: %s',string(datetime('now'))));
    writetable(best_exp_tab, fname,'WriteMode','append',"WriteVariableNames",true);

    fname = fullfile([datadir,pname,'_',ptype,'_',solveargs{1}{2},'_all','.txt']);
    fid=fopen(fname,'w+');
    fprintf(fid,'%s \n',sprintf('%% file: %s',mfilename));
    fprintf(fid,'%s \n',sprintf('%% date: %s',string(datetime('now'))));
    fclose(fid);
    writetable(all_exp_tab, fname, 'WriteMode','append',"WriteVariableNames",true);

    fname = fullfile([datadir,pname,'_',ptype,'_',solveargs{1}{2},'_trace']);
    save(fname,'runtrace');

end

fprintf('run end: %s\n',string(datetime('now')));
diary off;
