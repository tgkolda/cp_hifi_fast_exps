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

% Plot the original solution
ps1alt = @(x,y) plot(((size(Xorig,1)-1)./(length(x)-1)).*(x-1)+1,y,'.','Color',"#377EB8",'LineWidth',1,'MarkerSize',10);
ps2alt = @(x,y) plot(((size(Xorig,2)-1)./(length(x)-1)).*(x-1)+1,y,'.','Color',"#76C175",'LineWidth',1,'MarkerSize',10);
viz(Mtrue,'Figure',1,vizopts{:},'PlotCommands',{ps1alt,ps2alt,ps3},'Title', 'Full Data Solution');
image_resize(cpw,cph);

%% Create sampled tensors
nsamples = 2500./[1 2.5 5];
rng(0); 
seeds = (randi(1e9,length(nsamples),1));
for i = 1:length(nsamples)
    rng(seeds(i));
    X{i} = tensor_unaligned(Xorig,nsamples(i));
end

%% General Setup
nruns = 5;

%% Setup for CP-HIFI

rkhsinfo{1}.inf = true;
rkhsinfo{1}.kfunc = kernfunc_gaussian(15); 
rkhsinfo{1}.lambda = 0.1;
rkhsinfo{2}.inf = true;
rkhsinfo{2}.kfunc = kernfunc_gaussian(8); 
rkhsinfo{2}.lambda = 0.1;
rkhsinfo{3}.inf = false;

%% Run Loop for CP-HIFI
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
    viz(Mhifi{i}, 'Figure', 2+(i-1), vizopts{:}, 'title', sprintf('CP-HIFI with %d samples',nsamples(i)))
    image_resize(cpw,cph);
end

%% Compare to CP-WOPT without lower bound
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
        [Mtmp,~,info] = cp_wopt(Xfull/scalefactor,W,3);
        if j == 1 || (info.f < bestf)
            bestf = info.f;
            Mopt{i} = scalefactor*Mtmp;
        end
    end
    viz(Mopt{i}, 'Figure', 5+(i-1), vizopts{:}, 'title', sprintf('CP-WOPT with %d samples',nsamples(i)))
    image_resize(cpw,cph);
end

%% Compare to CP-WOPT with lower bound
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
        [Mtmp,~,info] = cp_wopt(Xfull/scalefactor,W,3, 'lower', 0);
        if j == 1 || (info.f < bestf)
            bestf = info.f;
            Moptlb{i} = scalefactor*Mtmp;
        end
    end
    viz(Moptlb{i}, 'Figure', 9+(i-1), vizopts{:}, 'title', sprintf('CP-WOPT with lower bound with %d samples',nsamples(i)))
    image_resize(cpw,cph);
end


%%
diary off