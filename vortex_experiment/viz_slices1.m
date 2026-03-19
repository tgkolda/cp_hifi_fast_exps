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