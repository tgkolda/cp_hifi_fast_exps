function [best_runs, all_runs, traces] = run_aligned(X,R,hifiinfo,genargs,solveargs,nruns,varargin)
%RUN_ALIGNED Runs multiple CP-ALS-HIFI solves on an aligned tensor
%
%   BEST_RUNS = RUN_ALIGNED(X, R, HIFIINFO, GENARGS, SOLVEARGS,
%   NRUNS) creates an aligned tensor from the full tensor X.
%   It then runs CP-ALS-HIFI NRUNS times for each solver
%   specified in SOLVEARGS, using the options in GENARGS for all runs. The
%   rank of the CP decomposition is R and the HIFI information is given in
%   HIFIINFO. The function returns BEST_RUNS, a table containing the best
%   run for each solver in terms of relative error.
%
%   BEST_RUNS = RUN_ALIGNED(X, R, HIFIINFO, GENARGS, SOLVEARGS,
%   NRUNS,SEEDS) specifies the seeds
%
%   BEST_RUNS = RUN_ALIGNED(X, R, HIFIINFO, GENARGS, SOLVEARGS,
%   NRUNS,SEEDS,DOSAVE) Saves the best hifi tensor
%
%   BEST_RUNS = RUN_ALIGNED(X, R, HIFIINFO, GENARGS, SOLVEARGS,
%   NRUNS,SEEDS,DOSAVE,SAVEPTH) Specifies the path to save the best tensors
%
%   [BEST_RUNS, ALL_RUNS] = RUN_ALIGNED(...) also returns ALL_RUNS, a table
%   containing the results of all runs.
%
%   [BEST_RUNS, ALL_RUNS, TRACES] = RUN_ALIGNED(...) also returns TRACES, a 
%   cell array with iteration information

% Johannes Brust          johannesbrust@yahoo.com
% School of Mathematical and Statistical Sciences
% Arizona State University

all_runs = table('Size',[length(solveargs)*nruns 7], ...
    'VariableTypes',{'string','double','double','double','double','double','double'}, ...
    'VariableNames',{'Solver','RunId','Seed','Time','RelErr','Fit','Iters'});

if 6<nargin
    seeds = varargin{1};
else
    seeds = 0:(nruns-1);
end
if 7<nargin
    dosave = varargin{2};
else
    dosave = 0;
end
if 8<nargin
    savepth = varargin{3};
else
    savepth = '.';
end
if dosave
    MS = cell(length(solveargs),2);
end

traces = cell(nruns,3);

for i = 1:length(solveargs)
    relerr_     = 1e3;
    t_elapsed_  = 0;
    id_         = 1;
    if dosave
        MS{i,1} = solveargs{i}{2};
    end
    for runid = 1:nruns
        fprintf('\n Run #%d for solver (using rng seed %d): %s\n', runid, seeds(runid), solveargs{i}{2});
        rng(seeds(runid),'twister');
        tic
        [M,info] = cp_als_hifi(X, R, hifiinfo, genargs{:}, solveargs{i}{:});
        t_elapsed = toc;
        relerr = norm(full(M)-X)/norm(X);
        all_runs((i-1)*nruns + runid,:) = ...
    		{solveargs{i}{2}, runid, seeds(runid), t_elapsed, relerr, info.final_fit, info.iters};

        [traces{runid,:}] = deal(runid,seeds(runid),info);

        if relerr < relerr_
            relerr_ = relerr;
            t_elapsed_ = t_elapsed;
            id_ = runid;
            if dosave
                MS{i,2} = M;
            end
        end

    end
    fprintf('        Solver=%s (best run at id=%i, seed=%i)\n', solveargs{i}{2},id_,seeds(id_));
    fprintf('Best rel error=%5.4f\n', relerr_);
    fprintf('     Best time=%5.4f\n', t_elapsed_);

end

% Find best run per solver in terms of RelErr
best_runs = extract_best_runs(all_runs, 'RelErr');
if dosave
    save(savepth,'MS');
end
