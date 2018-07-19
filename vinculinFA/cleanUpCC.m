function [CC allFeatures] = cleanUpCC(movieData, CC, allFeatures, ...
    tFirst, tLast, pFirst, bandWidth, alpha)

nFrames = movieData.nImages(1);
pixelSize = movieData.pixelSize_nm;
imSize = movieData.imSize;

bandWidth = bandWidth / pixelSize;
nCC = numel(CC);

distToEdgePath = movieData.distanceTransform.directory;
distToEdgeFiles = dir([distToEdgePath filesep '*.mat']);

% Compute the first and last frame of each CC
tFirstCC = cellfun(@(trackIdx) min(tFirst(trackIdx)), CC); % first frame of CC
tLastCC = cellfun(@(trackIdx) max(tLast(trackIdx)), CC);   % last frame of CC

N = zeros(nCC,1);
K = zeros(nCC,1);

% angular accuracy
p = .5;
dt = p * pi / 2;

for iFrame = 1:nFrames
    % Load distance transform
    load(fullfile(distToEdgePath, distToEdgeFiles(iFrame).name));

    isInFrame = find(iFrame >= tFirstCC & iFrame <= tLastCC);
    
    pFirstCC = cellfun(@(trackIdx) pFirst(trackIdx) + iFrame - tFirst(trackIdx), ...
        CC(isInFrame), 'UniformOutput', false);
    
    isTrackInFrame = cellfun(@(trackIdx) iFrame >= tFirst(trackIdx) & ...
        iFrame <= tLast(trackIdx), CC(isInFrame), 'UniformOutput', false);
    
    % Gather every feature of each track in the CC between at iFrame
    allFeaturesCC = cellfun(@(aa,bb) arrayfun(@(a,b) ...
        allFeatures(a:a+b-1, [1 2 4 6]), aa, bb, 'UniformOutput', false), ...
        pFirstCC, isTrackInFrame, 'UniformOutput', false);
    
    allFeaturesCC = cellfun(@(c) vertcat(c{:}), allFeaturesCC, ...
        'UniformOutput', false);
    
    modelsCC = getSegmentModels(allFeaturesCC);
    
    x = modelsCC(:,1);
    y = modelsCC(:,2);
    l = modelsCC(:,3);
    t = modelsCC(:,4);
    
    % Check that each segment is in the band...
    ind = sub2ind(size(distToEdge), round(y), round(x));
    isInBand = distToEdge(ind) < bandWidth;
    
    % ... and crop arrays accordingly.
    isInFrame = isInFrame(isInBand);
    x = x(isInBand);
    y = y(isInBand);
    l = l(isInBand);
    t = t(isInBand);

    % Compute extremity points
    ct = cos(t);
    st = sin(t);
    x1 = x + ct .* l/2;
    x2 = x - ct .* l/2;
    y1 = y + st .* l/2;
    y2 = y - st .* l/2;
    coords = [x1 y1 x2 y2];
    
    % Compute the unit vector along the pair
    dx = x2 - x1;
    dy = y2 - y1;
    norm = sqrt(dx.^2 + dy.^2);
    dx = dx ./ norm;
    dy = dy ./ norm;

    % Integer-valued extremity points
    iCoords = round(coords);

    % Compute Bresenham line between the 2 extremity points
    iCoords = num2cell(iCoords,1);
    pixels = arrayfun(@(ix1,iy1,ix2,iy2) bresenhamMEX([ix1,iy1], [ix2,iy2]), ...
        iCoords{:}, 'UniformOutput', false);

    % Remove pixels outside image boundaries
    pixels = cellfun(@(p) p(all(p > 0,2) & p(:,1) <= imSize(2) & ...
        p(:,2) <= imSize(1), :), pixels, 'UniformOutput', false);
    
    % Number of pixels in each segment
    nPixels = cellfun(@(p) size(p,1), pixels);
    
    % Concatenate every pixels
    ppLast = cumsum(nPixels);
    ppFirst = ppLast - nPixels + 1;
    pixels = vertcat(pixels{:});

    % Note: if numel(pixels) == 0, there will be an 'Index exceeds matrix
    % dimensions' error. To prevent that, we force pixels to be a Nx2
    % matrix.
    pixels = reshape(pixels, size(pixels,1), 2);
    ind = sub2ind(size(distToEdge), pixels(:,2), pixels(:,1));

    % Get cell edge orientation along the Bresenham lines
    [dU dV] = gradient(distToEdge);
    dU = dU(ind);
    dV = dV(ind);
    norm = sqrt(dU.^2 + dV.^2);
    
    % Some pixels can be outside the cell's footprint and therefore, the norm
    % of the gradient can be 0.
    isnnz = norm ~= 0;
    dU = dU(isnnz) ./ norm(isnnz);
    dV = dV(isnnz) ./ norm(isnnz);
    
    nValidPixels = arrayfun(@(a,b) nnz(isnnz(a:a+b-1)), ppFirst, nPixels);
    ppLast = cumsum(nValidPixels);
    ppFirst = ppLast - nValidPixels + 1;
    
    % Expand dx and dy
    ind = arrayfun(@(a,b) repmat(b,a,1), nValidPixels, (1:numel(dx))', ...
        'UniformOutput', false);
    ind = vertcat(ind{:});
     
    % Note: same rational than above.
    ind = reshape(ind,size(ind,1),1);
    dx = dx(ind);
    dy = dy(ind);
 
    % Compute the angle between the segment model and each gradient vectors
    theta = acos(abs(dx .* dU + dy .* dV));
     
    isAlongFlow = theta <= dt;
     
    N(isInFrame) = N(isInFrame) + nValidPixels;
    K(isInFrame) = K(isInFrame) + ...
        arrayfun(@(a,b) nnz(isAlongFlow(a:a+b-1)), ppFirst, nValidPixels);
end

% Binomial test: for each CC, we define the segment model. Among N pixels
% on each segment, we calculate how many K pixels were oriented in the same
% direction of the segment. If this number K is too low (i.e. the segment
% is basically oriented along the cell edge), we invalide the CC.
isInBand = N ~= 0;
cutoffs = icdf('bino', 1-alpha, (1:max(N))', p);

isValid = true(nCC,1);
isValid(isInBand) = K(isInBand) >= cutoffs(N(isInBand));

% Split every invalid CC
newCC = cellfun(@(trackIdx) arrayfun(@(t) {t}, trackIdx'), CC(~isValid), ...
    'UniformOutput', false);
% For each new CC, reset the length of detection to diffraction limit
newCC = vertcat(newCC{:});

for iCC = 1:numel(newCC)
    t = newCC{iCC};
    ind = pFirst(t):pFirst(t)+tLast(t)-tFirst(t);
    allFeatures(ind,4) = allFeatures(ind,5);
end

CC = [CC(isValid); newCC];


