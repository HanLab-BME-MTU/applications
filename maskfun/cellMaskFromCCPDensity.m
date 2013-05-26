function mask = cellMaskFromCCPDensity(frame, connect, showHist, modeRatio)

g = filterGauss2D(frame, 5);
[ny, nx] = size(frame);

% remove outliers
v = frame(:);
pct = prctile(v, [1 99]);
v(v<pct(1) | pct(2)<v) = [];
[f,xi] = ksdensity(v, 'npoints', 100, 'support', 'positive');

% local max/min
lmax = locmax1d(f, 3);
lmin = locmin1d(f, 3);

% max value
% hmax = find(f==max(f), 1, 'first');
dxi = xi(2)-xi(1);

% identify min after first mode
if ~isempty(lmin)
    idx = find(lmin>lmax(1), 1, 'first');
    if ~isempty(idx) && sum(f(1:lmin(idx(1))))*dxi < modeRatio
        min0 = lmin(idx);
        T = xi(min0);
        mask = g>T;
    else
        mask = ones(ny,nx);
    end
else
    mask = ones(ny,nx);
end

borderIdx = [1:ny (nx-1)*ny+(1:ny) ny+1:ny:(nx-2)*ny+1 2*ny:ny:(nx-1)*ny];

% retain largest connected component
if connect
    CC = bwconncomp(mask, 8);
    compsize = cellfun(@(i) numel(i), CC.PixelIdxList);
    mask = zeros(ny,nx);
    mask(CC.PixelIdxList{compsize==max(compsize)}) = 1;
end

% fill holes (retain largest boundary)
B = bwboundaries(mask);
nb = cellfun(@(i) numel(i), B);
B = B{nb==max(nb)};
boundary = zeros(ny,nx);
boundary(sub2ind([ny nx], B(:,1), B(:,2))) = 1;
% boundary(borderIdx) = 1;
CC = bwconncomp(1-boundary, 4);

% mask indexes
labels = double(labelmatrix(CC));
idx = unique(mask.*labels);
boundaryLabel = unique(labels(boundary==1));
% average intensities for each component
compLabels = setdiff(idx, boundaryLabel);
nLabels = numel(compLabels);
if nLabels > 1
    meanIntensity = zeros(1,nLabels);
    for k = 1:nLabels
        meanIntensity(k) = mean(frame(labels==compLabels(k)));
    end
    bgLabel = compLabels(find(meanIntensity==min(meanIntensity),1,'first'));
    mask = ismember(labels, setdiff(idx, bgLabel)); % 3rd version
else
    mask = ismember(labels, compLabels) | boundary;
end
% mask = boundary | ismember(labels, idx(2:end)); % 2nd version
% mask = boundary | labels==idx(2); % 1st version


if showHist
    dx = xi(2)-xi(1);
    ni = hist(v, xi);
    ni = ni/sum(ni)/dx;
    
    figure;
    hold on;
    plot(xi, ni, 'k.-', 'LineWidth', 3, 'MarkerSize', 20);
    plot(xi, f, 'r-', 'LineWidth', 1);
    set(gca, 'LineWidth', 2, 'FontSize', 18);
end