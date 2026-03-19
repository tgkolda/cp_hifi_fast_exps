function best_runs = extract_best_runs(all_runs, colname)
% EXTRACT_BEST_RUNS Extract the best run per solver from all runs table.
%
%   BEST_RUNS = EXTRACT_BEST_RUNS(ALL_RUNS) takes a table ALL_RUNS containing
%   results from multiple runs of different solvers and extracts the best run
%   for each solver based on the 'RelErr' metric.
%
%   BEST_RUNS = EXTRACT_BEST_RUNS(ALL_RUNS, COLNAME) takes a table ALL_RUNS
%   containing results from multiple runs of different solvers and extracts
%   the best run for each solver based on the metric specified by COLNAME.

if nargin < 2
	colname = 'RelErr';
end
G = findgroups(all_runs.Solver);
rowIdx = splitapply(@(r, idx) idx(find(r == min(r), 1)), all_runs.(colname), (1:height(all_runs))', G);
best_runs = all_runs(sort(rowIdx), :);
%best_runs = all_runs(rowIdx, :);
