% Example:
% results = analyze_collagen_fibers('your_sem_image.png', 5.0); % 5 nm/pixel
% results.DiameterSummary, results.CurvatureSummary, results.Persistence
function results = analyze_collagen_fibers(imPath, pixelSize_nm)
% Analyze collagen fibers in an SEM image.
% Inputs:
%   imPath        : path to image (8-bit/16-bit ok)
%   pixelSize_nm  : pixel size in nanometers (e.g., 5.0 nm/pixel)
% Outputs (struct):
%   DiameterPerPoint (table) : diameter at each skeleton point
%   CurvaturePerPoint (table): curvature at each arc-length sample
%   FiberTable (table)       : per-fiber summaries (length, mean diameter, mean curvature)
%   DiameterSummary          : mean/median/IQR (nm)
%   CurvatureSummary         : mean/median/IQR (1/µm)
%   Persistence              : estimated persistence length Lp (µm), fit diagnostics

%% --- 0. Read & normalize
I0 = im2double(imread(imPath));
if size(I0,3)>1, I0 = rgb2gray(I0); end

%% --- 1. Enhance ridges (SEM fibers are bright-on-dark or dark-on-bright?auto decide)
I = imgaussfilt(I0, 1.0);
if mean(I(:)) < 0.5
    % bright fibers on dark
    Ienh = I;
else
    % dark fibers on bright ? invert
    Ienh = imcomplement(I);
end

% Hessian/Frangi-like ridge enhancement (lightweight)
Ih = ridge_enhance(Ienh, 1:3);    % multi-scale (? = 1..3 px)
Ih = mat2gray(Ih);

%% --- 2. Segment & clean
T = adaptthresh(Ih, 0.45, 'NeighborhoodSize', 41, 'Statistic','gaussian');
BW = imbinarize(Ih, T);
BW = bwareaopen(BW, 50);                       % remove tiny specks
BW = imclose(BW, strel('disk', 1));            % fill small gaps
BW = imfill(BW, 'holes');

%% --- 3. Skeleton & local diameter via EDT
Skel = bwskel(BW, 'MinBranchLength', 10);      % 1-px centerlines
D = bwdist(~BW);                                % Euclidean distance
[ry, rx] = find(Skel);
r_px = D(Skel);                                 % radius at skeleton pixels (px)
d_nm = 2 * r_px * pixelSize_nm;                 % local diameter (nm)

% Connected components on skeleton ? individual fibers
CC = bwconncomp(Skel, 8);
label = labelmatrix(CC);

%% --- 4. Trace ordered centerlines for each fiber
fibers = cell(CC.NumObjects,1);
for k = 1:CC.NumObjects
    idx = CC.PixelIdxList{k};
    [y, x] = ind2sub(size(Skel), idx);
    [x_ord, y_ord] = order_skeleton_path(x, y);
    fibers{k} = [x_ord(:), y_ord(:)];          % [x, y] in pixels
end

%% --- 5. Per-fiber metrics: length, curvature, diameter distribution
px2um = pixelSize_nm / 1000;                    % µm per pixel
allCurv = []; allDiam = []; fpID = []; sCoord = []; tangents = {};
fiberLen_um = zeros(CC.NumObjects,1);
meanDiam_nm = zeros(CC.NumObjects,1);
meanCurv_invum = zeros(CC.NumObjects,1);

for k = 1:CC.NumObjects
    P = fibers{k};
    if size(P,1) < 10, continue; end

    % local diameter along this fiber (use nearest skeleton pixels)
    linIdx = sub2ind(size(Skel), round(P(:,2)), round(P(:,1)));
    diam_k = 2 * D(linIdx) * pixelSize_nm;     % nm

    % arc length in µm
    ds = [0; hypot(diff(P(:,1)), diff(P(:,2)))] * px2um;
    s = cumsum(ds);

    % tangent & curvature (finite difference on smoothed centerline)
    Ps = smoothdata(P, 1, 'rlowess', 9);       % light smoothing
    dP  = diff(Ps,1,1);                         % first differences
    dsP = hypot(dP(:,1), dP(:,2)) * px2um + eps;
    tx = dP(:,1) ./ (dsP/px2um); ty = dP(:,2) ./ (dsP/px2um);
    t  = [tx ty]; t = t ./ vecnorm(t,2,2);

    dt  = diff(t,1,1);
    ds2 = (s(3:end) - s(2:end-1));             % centered spacing (µm)
    kappa = vecnorm(dt,2,2) ./ max(ds2, eps);  % 1/µm, curvature magnitude
    % align lengths
    s_mid = s(2:end-1);
    diam_mid = diam_k(2:end-1);

    % summaries
    fiberLen_um(k)      = s(end);
    meanDiam_nm(k)      = median(diam_k,'omitnan');
    meanCurv_invum(k)   = median(kappa,'omitnan');

    % store for global tables
    allCurv = [allCurv; kappa];
    allDiam = [allDiam; diam_mid];
    fpID    = [fpID; k*ones(numel(kappa),1)];
    sCoord  = [sCoord; s_mid];
    tangents{k} = t; % for persistence
end

%% --- 6. Persistence length via tangent?tangent correlation (2D WLC)
% C(?s) = < t(s)·t(s+?s) > ? exp(-?s / Lp)  for 2D projection
% Build correlation over all fibers, binned by ?s
[dsBins, Cmean] = tangent_correlation(tangents, px2um, 1.0); % 1 µm bins
valid = ~isnan(Cmean) & dsBins>0 & Cmean>0 & Cmean<1;
x = dsBins(valid)'; y = Cmean(valid)';

% robust log-linear fit: y ? exp(-x/Lp)
beta = robustfit(x, log(y), 'huber');
Lp_um = -1 / beta(2);
yfit = exp(beta(1) + beta(2)*x);

%% --- 7. Pack outputs
DiameterPerPoint = table(fpID, sCoord, allDiam, ...
    'VariableNames', {'FiberID','s_um','Diameter_nm'});
CurvaturePerPoint = table(fpID, sCoord, allCurv, ...
    'VariableNames', {'FiberID','s_um','Curvature_per_um'});

FiberTable = table((1:CC.NumObjects)', fiberLen_um, meanDiam_nm, meanCurv_invum, ...
    'VariableNames', {'FiberID','Length_um','MedianDiameter_nm','MedianCurvature_per_um'});

results.DiameterPerPoint   = DiameterPerPoint;
results.CurvaturePerPoint  = CurvaturePerPoint;
results.FiberTable         = FiberTable;
results.DiameterSummary    = stats_summary(DiameterPerPoint.Diameter_nm);
results.CurvatureSummary   = stats_summary(CurvaturePerPoint.Curvature_per_um);
results.Persistence.Lp_um  = Lp_um;
results.Persistence.Bins_um = x;
results.Persistence.C_mean  = y;
results.Persistence.C_fit   = yfit;

% quick plots
figure; histogram(DiameterPerPoint.Diameter_nm, 60); xlabel('Diameter (nm)'); ylabel('Count'); title('Diameter distribution');
figure; histogram(CurvaturePerPoint.Curvature_per_um, 60); xlabel('Curvature (1/\mum)'); ylabel('Count'); title('Curvature distribution');
figure; plot(x, y, 'o'); hold on; plot(x, yfit, '-', 'LineWidth',2); 
xlabel('\Delta s (\mum)'); ylabel('<t(0)\cdot t(\Delta s)>'); title(sprintf('Persistence length L_p = %.2f \\mum', Lp_um));
grid on;

end

%% ---------- helpers ----------
function Ih = ridge_enhance(I, sigmas)
% Simple multi-scale ridge enhancement using Hessian eigenvalues.
I = im2double(I);
Ih = zeros(size(I));
for s = sigmas
    G = imgaussfilt(I, s);
    [L1, L2] = hessian_eigs(G, s);
    R = max(0, -(L2));  % prefer line-like (one large neg eigenvalue)
    Ih = max(Ih, R);
end
end

function [L1, L2] = hessian_eigs(I, s)
% Hessian eigenvalues at scale s (pixels). Sign convention for bright ridges.
[Ix, Iy] = gradient(I);
[Ixx, Ixy] = gradient(Ix);
[~,  Iyy]  = gradient(Iy);
% scale-normalized Hessian
Ixx = (s^2)*Ixx; Ixy = (s^2)*Ixy; Iyy = (s^2)*Iyy;
% eigenvalues of 2x2 symmetric matrix
tmp = sqrt( (Ixx - Iyy).^2 + 4*Ixy.^2 );
L1 = 0.5*(Ixx + Iyy - tmp);
L2 = 0.5*(Ixx + Iyy + tmp);
end

function [xs, ys] = order_skeleton_path(x, y)
% Order 8-connected skeleton pixels into a single path (handles branches by taking the longest geodesic).
P = [x(:) y(:)];
% build adjacency on pixel grid
[Xu, Yu] = ndgrid(min(x):max(x), min(y):max(y));
gridIdx = containers.Map('KeyType','char','ValueType','int32');
for i=1:size(P,1)
    key = sprintf('%d_%d', P(i,1), P(i,2));
    gridIdx(key) = i;
end
% compute degree
deg = zeros(size(P,1),1);
nbrs = cell(size(P,1),1);
K = [-1 0 1];
for i=1:size(P,1)
    xi=P(i,1); yi=P(i,2);
    nb = [];
    for dx=K, for dy=K
        if dx==0 && dy==0, continue; end
        key = sprintf('%d_%d', xi+dx, yi+dy);
        if isKey(gridIdx, key), nb = [nb gridIdx(key)]; end
    end, end
    nbrs{i} = nb; deg(i) = numel(nb);
end
ends = find(deg==1);
if isempty(ends)
    % loop: pick arbitrary start and unwrap by DFS
    start = 1; [path] = walk(nbrs, start);
else
    % choose furthest pair (approximate longest path)
    [path1] = walk(nbrs, ends(1));
    start = path1(end);
    [path] = walk(nbrs, start);
end
xs = P(path,1); ys = P(path,2);
end

function path = walk(nbrs, start)
N = numel(nbrs);
visited = false(N,1);
path = zeros(N,1); k=0;
stack = [start 0]; % [node, parent]
while ~isempty(stack)
    node = stack(end,1); parent = stack(end,2);
    stack(end,:) = [];
    if visited(node), continue; end
    visited(node)=true; k=k+1; path(k)=node;
    for nb = nbrs{node}
        if nb==parent, continue; end
        if ~visited(nb), stack = [stack; nb node]; end %#ok<AGROW>
    end
end
path = path(1:k);
end

function [dsBins, Cmean] = tangent_correlation(tangents, px2um, bin_um)
% Build tangent?tangent correlation vs contour separation ?s.
pairs = {};
maxS = 0;
for k=1:numel(tangents)
    t = tangents{k};
    if size(t,1)<5, continue; end
    pairs{end+1} = t; %#ok<AGROW>
    maxS = max(maxS, (size(t,1)-1)*px2um);
end
edges = 0:bin_um:(maxS+bin_um);
numBins = numel(edges)-1;
Csum = zeros(1,numBins); Ccnt = zeros(1,numBins);
for k=1:numel(pairs)
    t = pairs{k};
    M = size(t,1);
    for lag = 1:min(300, M-1) % up to ~300 steps to keep it fast
        ds = lag*px2um;
        b = discretize(ds, edges);
        if isnan(b), continue; end
        dots = sum(t(1:end-lag,:).*t(1+lag:end,:), 2);
        Csum(b) = Csum(b) + sum(dots);
        Ccnt(b) = Ccnt(b) + numel(dots);
    end
end
Cmean = Csum ./ max(Ccnt,1);
Cmean(Ccnt==0) = NaN;
dsBins = (edges(1:end-1)+edges(2:end))/2;
end

function S = stats_summary(x)
S.mean   = mean(x,'omitnan');
S.median = median(x,'omitnan');
S.IQR    = iqr(x);
S.std    = std(x,'omitnan');
end