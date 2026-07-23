%% Run Amino Acids Experiments (with Missing Data)

%% Setup paths
if ~exist('tensor.m','file')
    fprintf('Adding tensor toolbox to path...\n');
    addpath('../../../tensor_toolbox/');
end
if ~exist('tensor.m','file')
    error('Tensor Toolbox not found.');
end
if ~exist('cp_als_hifi.m','file')
    fprintf('Adding CP-HiFi code to path...\n');
    addpath('../../cp_hifi_code/');
end
if ~exist('cp_als_hifi.m','file')
    error('CP-HiFi code not found.');
end

%% Load the data
load(fullfile(getfield(what('tensor_toolbox'),'path'),'doc','aminoacids.mat'))
clear aminoacids_url aminoacids_cite
X = permute(X,[2 3 1]);
xrange = [min(X(:)), max(X(:))];
Xorig = X;
clear X

%% General experiment parameters (defined early so the diary header can report them)
nsamples = 2500./[0.5 1 2.5 5];
nruns = 5;

%% Setup results directory and start diary
resultsdir = '../results/';
if ~exist(resultsdir,'dir'); mkdir(resultsdir); end
pname = 'amino';
rr = matlabRelease;
ss = strsplit(char(datetime('now')));
dn = [ss{1},'-',strrep(ss{2},':','_')];
diary([resultsdir,pname,'_log_',dn,'.txt']);

fprintf('\n*************** %s experiment **********          \n',pname);
fprintf('*           Software: Matlab %s                     \n',rr.Release);
fprintf('*             System: %s                            \n',computer);
fprintf('*               File: %s                            \n',mfilename);
fprintf('*                                                   \n');
fprintf('*         nruns(ALS): %s                            \n',num2str(nruns));
fprintf('*           nsamples: %s                            \n',num2str(nsamples));
fprintf('*                                                   \n');
fprintf('*           datetime: %s                            \n',string(datetime('now')));
fprintf('*************************************************** \n\n');

%% Option: save figures to disk
% Set savefigs = false to skip saving. Figures are written as PNGs
% into figdir (default: the current amino_acids directory), one per solution.
savefigs = true;
figdir = './';
savefig_fn = @(nm) exportgraphics(gcf, [figdir, nm, '.png'], ...
    'Resolution', 300);


%% Setup for MATLAB plotting
% The plot function in the third mode will be different in each case
ps3 = @(x,y) bar(x,y,'FaceColor',"#984DA3");
ps1 = @(x,y) plot(((size(Xorig,1)-1)./(length(x)-1)).*(x-1)+1,y,'Color',"#377EB8",'LineWidth',2);
ps2 = @(x,y) plot(((size(Xorig,2)-1)./(length(x)-1)).*(x-1)+1,y,'-','Color',"#76C175",'LineWidth',2);
vizopts = {...
	'ModeTitles',{'Emission','Excitation','Sample'}, ...
    'FactorTitles','None', ...
	'XTicks',false(3,1),...
	'Normalize', @(T) arrange(T,1),...    
    'RelModeWidth', [0.5 0.3 0.2], ...
	'BaseFont',18, ...
    'TopSpace',0.15,...
	'BottomSpace',0.02,...
	'LeftSpace',0.02,...
    'Xlims',{[0 size(Xorig,1)+1],[0 size(Xorig,2)+1],[0 size(Xorig,3)+1]},...
	'PlotCommands',{ps1,ps2,ps3}
	};
cpw = 7;
cph = 5;

%% Get the "true" solution, just for comparisons later.
rng('default');
Mtrue = cp_als(Xorig,3);
REtrue = norm(full(Mtrue)-Xorig)/norm(Xorig);
fprintf('Rel error for True: %g\n', REtrue);

% Write out the true-solution factor matrices for plotting
write_factor_dat_files(Mtrue, [resultsdir,'amino_true']);

% Plot the original solution
ps1alt = @(x,y) plot(((size(Xorig,1)-1)./(length(x)-1)).*(x-1)+1,y,'.','Color',"#377EB8",'LineWidth',1,'MarkerSize',10);
ps2alt = @(x,y) plot(((size(Xorig,2)-1)./(length(x)-1)).*(x-1)+1,y,'.','Color',"#76C175",'LineWidth',1,'MarkerSize',10);
viz(Mtrue,'Figure',1,vizopts{:},'PlotCommands',{ps1alt,ps2alt,ps3},'Title', 'Full Data Solution');
image_resize(cpw,cph);
if savefigs, savefig_fn('amino_true'); end

%% Create sampled tensors
rng(42);
seeds = (randi(1e9,length(nsamples),1));
X = cell(length(nsamples),1);
for i = 1:length(nsamples)
    rng(seeds(i));
    X{i} = tensor_unaligned(Xorig,nsamples(i));
end

%% Setup for CP-HIFI

rkhsinfo{1}.inf = true;
rkhsinfo{1}.kfunc = kernfunc_gaussian(15); 
rkhsinfo{1}.lambda = 0.1;
rkhsinfo{2}.inf = true;
rkhsinfo{2}.kfunc = kernfunc_gaussian(8); 
rkhsinfo{2}.lambda = 0.1;
rkhsinfo{3}.inf = false;

%% Run Loop for CP-HIFI
Mhifi_orig = cell(length(nsamples),1);
Mhifi = cell(length(nsamples),1);
for i = 1:length(nsamples)
    fprintf('--- CP-HIFI with %d samples ---\n\n',nsamples(i))
    relerr1 = 1e4; % init to something huge
    for j = 1:nruns
        fprintf('Run %d\n',j)
        Mtmp = cp_als_hifi(X{i},3,rkhsinfo);
        tmprelerr = norm_masked_diff(Mtmp,X{i})/norm(X{i});
        if tmprelerr < relerr1
            relerr1 = tmprelerr;
            Mhifi{i} = Mtmp;
        end
        if relerr1 < .04
            break
        end
    end
    Mhifi_orig{i} = Mhifi{i};
    Mhifi{i} = resample_modes(Mhifi{i}, size(Xorig));
    viz(Mhifi{i}, 'Figure', 2+(i-1), vizopts{:}, 'title', sprintf('CP-HIFI with %d samples',nsamples(i)))
    image_resize(cpw,cph);
    if savefigs, savefig_fn(sprintf('amino_hifi_%d',nsamples(i))); end
    write_factor_dat_files(Mhifi{i}, sprintf('%samino_hifi_%d',resultsdir,nsamples(i)));
end

%% Compare to CP-WOPT without lower bound
Mopt_orig = cell(length(nsamples),1);
Mopt = cell(length(nsamples),1);
for i = 1:length(nsamples)
    fprintf('--- CP-WOPT with %d samples ---\n\n',nsamples(i))
    W = tenzeros(size(X{i}));
    W(X{i}.subs) = 1;
    Xfull = tenzeros(size(X{i}));
    Xfull(X{i}.subs) = X{i}.vals;
    scalefactor = norm(X{i});
    for j = 1:nruns
        fprintf('Run %d\n',j)
        [Mtmp,~,info] = cp_wopt(Xfull/scalefactor,W,3,'verbosity',500);
        if j == 1 || (info.f < bestf)
            bestf = info.f;
            Mopt{i} = scalefactor*Mtmp;
        end
    end
    % Convert to ktensor_hifi (preserves the x-value mapping)
    Mopt_orig{i} = Mopt{i};
    Mopt{i} = ktensor_hifi(Mopt{i}.lambda, Mopt{i}.U, X{i}.xvals, {[],[],[]}, {[],[],[]});
    viz(Mopt{i}, 'Figure', 5+(i-1), vizopts{:}, 'title', sprintf('CP-WOPT with %d samples',nsamples(i)))
    image_resize(cpw,cph);
    if savefigs, savefig_fn(sprintf('amino_wopt_%d',nsamples(i))); end
    write_factor_dat_files(Mopt{i}, sprintf('%samino_wopt_%d',resultsdir,nsamples(i)));
end

%% Compare to CP-WOPT with lower bound
Moptlb_orig = cell(length(nsamples),1);
Moptlb = cell(length(nsamples),1);
for i = 1:length(nsamples)
    fprintf('--- CP-WOPT with lower bound and %d samples ---\n\n',nsamples(i))
    W = tenzeros(size(X{i}));
    W(X{i}.subs) = 1;
    Xfull = tenzeros(size(X{i}));
    Xfull(X{i}.subs) = X{i}.vals;
    scalefactor = norm(X{i});
    for j = 1:nruns
        fprintf('Run %d\n',j)
        [Mtmp,~,info] = cp_wopt(Xfull/scalefactor,W,3, 'lower', 0,'verbosity',500);
        if j == 1 || (info.f < bestf)
            bestf = info.f;
            Moptlb{i} = scalefactor*Mtmp;
        end
    end
    Moptlb_orig{i} = Moptlb{i};
    Moptlb{i} = ktensor_hifi(Moptlb{i}.lambda, Moptlb{i}.U, X{i}.xvals, {[],[],[]}, {[],[],[]});
    viz(Moptlb{i}, 'Figure', 9+(i-1), vizopts{:}, 'title', sprintf('CP-WOPT with lower bound with %d samples',nsamples(i)))
    image_resize(cpw,cph);
    if savefigs, savefig_fn(sprintf('amino_woptlb_%d',nsamples(i))); end
    write_factor_dat_files(Moptlb{i}, sprintf('%samino_woptlb_%d',resultsdir,nsamples(i)));
end


%%
diary off