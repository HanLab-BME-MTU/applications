function results = analyze_collagen_sem(imPath, pixelSize_nm)
% Centerline + diameter, curvature & persistence from an SEM fiber image.
% Requires Image Processing Toolbox (uses fibermetric, bwskel).
% pixelSize_nm: set [] if unknown; outputs default to pixel units.

if nargin < 2, pixelSize_nm = []; end
I0 = im2double(imread(imPath));
if size(I0,3)>1, I0 = rgb2gray(I0); end

% --- 1) Contrast + steerable ridge enhancement (multi-scale)
I = adapthisteq(I0,'clipLimit',0.01);
sigma = [1 2 3];                            % (px) adjust if fibers thinner/thicker
F = max(cat(3, fibermetric(I,sigma(1),'StructureSensitivity',1), ...
               fibermetric(I,sigma(2),'StructureSensitivity',1), ...
               fibermetric(I,sigma(3),'StructureSensitivity',1)), [], 3);

% --- 2) Segment and clean
T = graythresh(imgaussfilt(F,1));
BW = imgaussfilt(F,1) > T;
BW = bwareaopen(BW,80);
BW = imclose(BW, strel('disk',1));
BW = imfill(BW,'holes');

% --- 3) Skeleton and local diameter by distance transform
Skel = bwskel(BW, 'MinBranchLength',10);
D = bwdist(~BW);                        % radius in px on mask
[ry, rx] = find(Skel);
r_px = D(Skel); d_px = 2*r_px;

% Label individual fibers on skeleton
CC = bwconncomp(Skel,8);
L = labelmatrix(CC);

% --- 4) Trace ordered centerlines, curvature and per-fiber summaries
px2um = []; if ~isempty(pixelSize_nm), px2um = pixelSize_nm/1000; end
diam_rows = []; curv_rows = []; fiber_rows = [];

for k = 1:CC.NumObjects
    idx = CC.PixelIdxList{k};
    [y, x] = ind2sub(size(Skel), idx);
    [x_ord, y_ord] = order_path(x, y);           % below
    P = [x_ord(:) y_ord(:)];
    if size(P,1)<10, continue; end

    % arclength
    ds = [0; hypot(diff(P(:,1)), diff(P(:,2)))];
    s = cumsum(ds);

    % diameter along centerline
    linIdx = sub2ind(size(Skel), round(P(:,2)), round(P(:,1)));
    d_here = 2*D(linIdx);

    % curvature by sliding circle fit (win = 9 pts)
    [s_c, kappa] = curvature_polyline(P, 9);     % below

    % store rows
    diam_rows = [diam_rows; [k*ones(numel(s),1), s, d_here]];
    if ~isempty(kappa)
        keep = ~isnan(kappa);
        curv_rows = [curv_rows; [k*ones(nnz(keep),1), s_c(keep), kappa(keep)]];
        med_kappa = median(kappa,'omitnan');
    else
        med_kappa = NaN;
    end
    fiber_rows = [fiber_rows; k, s(end), median(d_here,'omitnan'), med_kappa];
end

Diameter = array2table(diam_rows, 'VariableNames',{'FiberID','s_px','Diameter_px'});
Curvature = array2table(curv_rows, 'VariableNames',{'FiberID','s_px','Curvature_per_px'});
FiberTbl = array2table(fiber_rows, 'VariableNames',{'FiberID','Length_px','MedianDiameter_px','MedianCurvature_per_px'});

% --- 5) Persistence length from tangent?tangent correlation
paths = extract_paths(Skel);                      % list of ordered [x,y] arrays
[ds_bins, Cmean] = tangent_corr(paths, 1, 400);   % below
valid = ds_bins>0 & Cmean>0 & Cmean<1 & ~isnan(Cmean);
x = ds_bins(valid); y = Cmean(valid);
mdl = fitlm(x, log(y), 'Intercept', true);
b = mdl.Coefficients.Estimate;                    % [a; b] where log y ? a + b x
Lp_px = -1/b(2);

% Optional unit conversion
if ~isempty(px2um)
    Diameter.Diameter_nm = Diameter.Diameter_px * pixelSize_nm;
    Curvature.Curvature_per_um = Curvature.Curvature_per_px / px2um;
    FiberTbl.Length_um = FiberTbl.Length_px * px2um;
    FiberTbl.MedianDiameter_nm = FiberTbl.MedianDiameter_px * pixelSize_nm;
    FiberTbl.MedianCurvature_per_um = FiberTbl.MedianCurvature_per_px / px2um;
    Lp_um = Lp_px * px2um;
else
    Lp_um = NaN;
end

% --- 6) Package results
results.DiameterPerPoint = Diameter;
results.CurvaturePerPoint = Curvature;
results.FiberTable = FiberTbl;
results.Persistence.Lp_px = Lp_px;
results.Persistence.Lp_um = Lp_um;
results.Persistence.fit_x_px = x;
results.Persistence.fit_y = y;

% quick plots
figure; histogram(Diameter.Diameter_px,60); xlabel('Diameter (px)'); ylabel('Count'); title('Diameter');
figure; histogram(Curvature.Curvature_per_px,60); xlabel('Curvature (1/px)'); ylabel('Count'); title('Curvature');
figure; plot(x, y, 'o'); hold on; plot(x, exp(b(1)+b(2)*x), '-', 'LineWidth',2);
xlabel('\Delta s (px)'); ylabel('<t\cdot t>'); title(sprintf('L_p \\approx %.1f px', Lp_px)); grid on;

end

% ===== Helpers =====
function [xs, ys] = order_path(x, y)
% Order 8-connected skeleton pixels into a single path (favoring the longest).
P = [x(:) y(:)];
H = max(y)-min(y)+1; W = max(x)-min(x)+1;
off = [min(x) min(y)]-1;                       % offset for indexing
A = false(H,W); A(sub2ind([H W], y-off(2), x-off(1))) = true;
deg = conv2(double(A), ones(3), 'same') - double(A);
ends = find(deg(:)==1);
B = bwlabel(A,8);
% build geodesic from an endpoint (or arbitrary if loop)
if ~isempty(ends)
    [ey, ex] = ind2sub([H W], ends(1));
else
    [ey, ex] = ind2sub([H W], find(A,1));
end
% simple walk
nbr = [-1 -1;-1 0;-1 1;0 -1;0 1;1 -1;1 0;1 1];
ys = []; xs = []; cur = [ey ex]; prev = [NaN NaN];
while A(cur(1),cur(2))
    ys(end+1)=cur(1); xs(end+1)=cur(2); %#ok<AGROW>
    A(cur(1),cur(2)) = false;
    moved = false;
    for i=1:8
        n = cur + nbr(i,:);
        if n(1)<1||n(1)>H||n(2)<1||n(2)>W, continue; end
        if A(n(1),n(2))
            cur = n; moved = true; break;
        end
    end
    if ~moved, break; end
end
xs = xs(:)+off(1); ys = ys(:)+off(2);
end

function [s, kappa] = curvature_polyline(P, win)
if nargin<2, win=9; end
if size(P,1) < win, s=[]; kappa=[]; return; end
ds = [0; hypot(diff(P(:,1)), diff(P(:,2)))];
s = cumsum(ds);
half = floor(win/2); kappa = nan(size(P,1),1);
for i = 1+half : size(P,1)-half
    pts = P(i-half:i+half, :);
    x = pts(:,1); y = pts(:,2);
    x = x-mean(x); y = y-mean(y);
    Suu = sum(x.^2); Svv = sum(y.^2); Suv = sum(x.*y);
    Suuu = sum(x.^3); Svvv = sum(y.^3);
    Suvv = sum(x.*y.^2); Svuu = sum(y.*x.^2);
    A = [Suu Suv; Suv Svv]; b = 0.5*[Suuu+Suvv; Svvv+Svuu];
    c = A\b; R = mean(sqrt((x-c(1)).^2 + (y-c(2)).^2));
    if R>0, kappa(i)=1/R; end
end
end

function paths = extract_paths(Skel)
CC = bwconncomp(Skel,8);
paths = cell(0,1);
for k=1:CC.NumObjects
    idx = CC.PixelIdxList{k}; [y,x] = ind2sub(size(Skel), idx);
    [xo, yo] = order_path(x,y);
    P = [xo yo]; if size(P,1)>=5, paths{end+1}=P; end %#ok<AGROW>
end
end

function [dsBins, Cmean] = tangent_corr(paths, step_px, maxLag)
if nargin<2, step_px=1; end
if nargin<3, maxLag=400; end
sums = zeros(maxLag,1); cnt = zeros(maxLag,1);
for i=1:numel(paths)
    P = paths{i};
    v = gradient(P); n = sqrt(sum(v.^2,2))+eps; t = v./n;
    L = size(t,1);
    for lag=1:min(maxLag-1,L-1)
        dots = sum(t(1:end-lag,:).*t(1+lag:end,:),2);
        sums(lag)=sums(lag)+sum(dots); cnt(lag)=cnt(lag)+numel(dots);
    end
end
Cmean = sums./max(cnt,1); Cmean(cnt==0)=NaN;
dsBins = (0:maxLag-1)'*step_px;
end