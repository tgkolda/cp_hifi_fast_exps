%%PLOT_VORTEX_SLICES_A Runs unaligned vortex experiment to
% reproduce the slices presented in the article
%
% See also CP_ALS_HIFI.

% Johannes Brust          johannesbrust@yahoo.com
% School of Mathematical and Statistical Sciences
% Arizona State University
%
% 01/09/26, J.B., reorganizing how to save data
% 01/18/26, J.B., reproduce vortex slices
% 01/29/26, J.B., Remove labelling
% 02/03/26, J.B., Use only exportgraphics

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
ptype   = 'unaligned';
datadir = '../paper_figs/';

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
hifiinfo{1}.rho = 1e-6; % 1e-6
hifiinfo{2}.rho = hifiinfo{1}.rho;
hifiinfo{3}.rho = hifiinfo{1}.rho;

%% Experiment setup
rng(0,'twister');
X        = Xorig;
R        = 25; % 25
nruns    = 1;
genargs  = {'printitn',10,'tol',1e-6,'maxiters',50,'trace',true};
%selidx   = 5;
nsamples = 50000;
rho      = 1e-6; % 0
xrange   = [min(X(:)),max(X(:))];
seeds    = [99971,99991];
hasdirect = 1;

% Print experiment infos and start diary
rr  = matlabRelease;
ss  = strsplit(char(datetime('now')));
dn  = [ss{1},'-',strrep(ss{2},':','_')];
%diary([datadir,pname,'_',ptype,'_','log','_',dn,'.txt']);

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


fprintf('\n Creating unaligned tensor with %d samples...\n', nsamples);
rng(0,'twister');
X       = tensor_unaligned(X, nsamples);
solvs   = {'pcg','direct_nonsym'};
solvs1  = {'pcg','direct nonsym'};
solvs   = {'pcg'};
solvs1  = {'pcg'};

MS      = cell(length(solvs),1);
times   = zeros(length(solvs),1);
errs    = zeros(length(solvs),1);

solveargs       = cell(2,1);
solveargs{1}    = {'solver','pcg','inmaxit', 75, 'intol', 1e-6};
solveargs{2}    = {'solver','direct_nonsym'};

%% Run CP_ALS_RKHS
% Timers and iteration counter 

normX   = norm(X);

for iter = 1:length(solvs)

    itt     = 0;
    itss    = 0;
    itm     = 1;
    ts      = 0;
    relerr  = 10;

    tcpals  = tic;

    seed = seeds(iter);
    rng(seed,'twister');
    solver  = solvs{iter};
    fprintf('\n cp-als run: %i \n',iter);
    fprintf('\n     solver: %s \n',solver);
    fprintf('     rng seed: %i \n',seed);

    %% alshifi with approximate solve
    ts_tmp = tic;
    [M_tmp,info] = cp_als_hifi(X,R,hifiinfo,genargs{:},solveargs{iter}{:});
    ts_tmp  = toc(ts_tmp);
    itss_tmp= info.iters;

    relerr_tmp = norm_masked_diff(M_tmp,X)/normX;
    if relerr_tmp < relerr
        M       = M_tmp;
        relerr  = relerr_tmp;
        ts      = ts_tmp;
        itss    = itss_tmp;
        itt     = seed;
        itm     = iter;
        if relerr_tmp < 1e-7 % 0.08 -- small error to not exit early
            break
        end
    end

    tcpals  = toc(tcpals);

    %% Plot the factors for the best result
    fprintf('Time (sec): %g \n', tcpals)
    fprintf('Time best rel error    : %g\n', ts)
    fprintf('Run best rel error     : %i\n', itm)
    fprintf('Best rel error for HIFI: %g\n', relerr)
    M_hifi = resample_modes(M,size(Xorig));

    % remove nans, which might be caused by "resample_modes" and
    % division by zero (some column norms may be zero)
    for i=1:length(M_hifi.u)
        U = M_hifi.u{i};
        U(isnan(U)) = 0;
        M_hifi.u{i} = U;
    end

    MS{iter} = M_hifi;
    times(iter) = ts;
    errs(iter) = relerr;

end

%% Plot slices

labels  = {'(a)','(c)','(d)','(b)'};
xanot   = 0.42;
xanot1  = xanot;
ylanot  = 0.34;
xlanot  = 0.14;
fntsz   = 8;

aligned = {'unaligned','unaligned'};

% plot original
figsize = [2.5,2.5];    % [4.5,4.5]
figure('PaperSize',figsize,'PaperUnits','inches');

% slice function
viz_slices1(Xorig,1,xrange,'nrow',1,'ncol',1);

% export figure
exportgraphics(gcf,[datadir,pname,'_original_a_exp.pdf'],'ContentType','vector','BackgroundColor','none');

% Plots
for i=1:length(errs)

    figsize = [2.5,2.5];    % [4.5,4.5]
    figure('PaperSize',figsize,'PaperUnits','inches');
    
    % slice function
    viz_slices1(MS{i},1,xrange,'nrow',1,'ncol',1);
    
    % export
    exportgraphics(gcf,[datadir,pname,'_slices_',solvs{i},'_a_exp.pdf'],'ContentType','vector','BackgroundColor','none');
            
end

% plot samples
figsize = [2.5,2.5];    % [4.5,4.5]
figure('PaperSize',figsize,'PaperUnits','inches');

% slice function
viz_slices1(X,1,xrange,'nrow',1,'ncol',1);

% annotations
xsz = size(X); %[256,256,7];
p = size(MS{1},3);
ttxt = ['$\textnormal{sampled } ',num2str(nsamples/prod(xsz)*100,1),...
    ' \% $'];
xtxt = ['$\textnormal{slice } ',num2str(p,3),'$'];
%title(ttxt,'interpreter','latex','FontSize',fntsz);
%xlabel(xtxt,'interpreter','latex','FontSize',fntsz,...
%    'Units', 'normalized', 'Position', [0.5, -0.2, 0]);

%la = annotation('textbox',[xlanot,ylanot,0.1,0.1],'String',labels{i+2},...
%        'Interpreter','Latex','EdgeColor','none','FontSize',fntsz-1);

hold on;
recpos = [1,1,xsz(2)-1,xsz(1)-1];
rectangle('Position', recpos, 'EdgeColor', 'k', 'LineWidth', 0.1);

%print([datadir,pname,'_samp_a'],'-dpdf','-fillpage');
exportgraphics(gcf,[datadir,pname,'_samp_a_exp.pdf'],'ContentType','vector','BackgroundColor','none');

%% local viz function
function xrange = viz_slices1(X,figid,xrange,varargin)
%% Set algorithm parameters from input or by using defaults
params = inputParser;
params.addParameter('nrows',5);
params.addParameter('ncols',2);
params.addParameter('cmap',colormap(jet));
params.parse(varargin{:});
nrows = params.Results.nrows;
ncols = params.Results.ncols;
cmap = params.Results.cmap;

if isa(X,'tensor')
    X = X.data;
elseif isa(X,'ktensor')
    X = full(X);
    X = X.data;
elseif isa(X,'tensor_unaligned')
    P = spones(X);
    P = double(full(P));
    X = double(full(X));
    X(P==0) = NaN;
    cmap = [1 1 1; cmap];
else
    error('Invalid tensor');
end

% if figid == 0
%     figure;
% else
%     figure(figid); clf;
% end

if nargin < 3
    xrange = [];
end

if isempty(xrange)
    xrange = [min(X(:)),max(X(:))];
end

p = size(X,3);
N = nrows * ncols;

idx = round(linspace(1,p,N));
idx(1) = 1;
idx(end) = p;

for i = 1:N
    j = idx(i);
    %subplot(nrows,ncols,i);
    imshow(X(:,:,j),xrange,'InitialMagnification','fit');
    colormap(cmap);
    %title(['Slice ' int2str(j)])
end

end