function movieData = getMoviePairTracks(movieData,iChannel,bandWidth,...
    minOverlap,timeGap,maxEuclidianDist,kSigma,alpha,probBinSize,...
    batchMode)

% iChannel              image channel where the tracks has been computed
%                       from.
%
% bandWidth             distance in nanometers away from the cell edge
%                       where pair of feature needs to be evaluate whether
%                       it follows the Actin flow.
%
% minOverlap:           is the minimal number of frame 2 tracks need to
%                       live together to be considered as potential pair of
%                       tracks. Default is 1.
%
% timeGap:              is the maximum number of frames 2 tracks can be
%                       separated from each other to be considered as
%                       potential pair of tracks. Default is 0.
%
% maxEuclidianDist:     maximum euclidian distance between the mean
%                       position of 2 tracks to be considered as potential
%                       pair of tracks.
%
% kSigma:               cutoff in number of standard deviations
%
% alpha:                quantile of PDF tail. Default is 0.05.
%
% probBinSize:          size in [0...1] of a bin. probBinSize^-1 gives the
%                       the number of bins a PDF is cut into. Default is
%                       1e-4.

initFunc = {...
    @initPairParams00,...
    @initPairParams01,...
    @(params1,params2) initPairParams01(params2,params1),...
    @initPairParams11};
    
modelSupportFunc = {...
    @anisoGaussian2DSupport,...
    @segment2DSupport};

modelFunc = {...
    @anisoGaussian2D,...
    @segment2D};

fitModelFunc = {...
    @fitAnisoGaussian2D,...
    @fitSegment2D};

%% Parse input parameters

if nargin < 3 || isempty(bandWidth)
    bandWidth = 1000;
end

if nargin < 4 || isempty(minOverlap)
    minOverlap = 1;
end

if nargin < 5 || isempty(timeGap)
    timeGap = 0;
end

if nargin < 8 || isempty(alpha)
    alpha = 0.05;
end

if nargin < 9 || isempty(probBinSize)
    probBinSize = 1e-4;
end

assert(checkMovieDistanceTransform(movieData));
assert(checkMovieParticleDetection(movieData));
assert(checkMovieParticleTracking(movieData));

movieData.pairTracks.status = 0;

movieData.pairTracks.directory = fullfile(movieData.analysisDirectory, 'pairTracks');

% Create output directory
if ~exist(movieData.pairTracks.directory, 'dir')
    mkdir(movieData.pairTracks.directory);
end

nFrames = movieData.nImages(1);
imSize = movieData.imSize;
ny = imSize(1);
nx = imSize(2);

%% Load feature and track infos

load(fullfile(movieData.particleDetection.directory, movieData.particleDetection.filename));
load(fullfile(movieData.particleTracking.directory, movieData.particleTracking.filename));

nTracks = size(tracksFinal,1); %#ok<NODEF>

% Check there is no split and merge. If split and merge is enabled, it
% would complexify the interpolation in gaps.
seqOfEvents = vertcat(tracksFinal(:).seqOfEvents);
assert(nnz(isnan(seqOfEvents(:,4))) == size(seqOfEvents,1));

SEL = getTrackSEL(tracksFinal);
tFirst = SEL(:,1);                   % first frame of the track
tLast = SEL(:,2);                    % last frame of the track
lifetime = SEL(:,3);                 % lifetime of the track

%% Linearize frame indices of all the tracks
allTrackFrameIdx = arrayfun(@(a,b) a:b, tFirst, tLast, 'UniformOutput', false);
allTrackFrameIdx = [allTrackFrameIdx{:}];

%% Linearize feature parameters of all tracks 

% allTrackParams: feature parameters of all tracks appended together
% x11 y11 A11 Sx11 Sy11 T11 C11
% x12 ... (2nd point of the 1st track)
% ...
% x1n ... (nth point of the 1st track)
% x21 ... (1st point of the 2nd track)
% ...
allTrackParams = nan(numel(allTrackFrameIdx),7);

tracksFeatIndx = [tracksFinal(:).tracksFeatIndxCG];

for iFrame = 1:nFrames
    ind = allTrackFrameIdx == iFrame & tracksFeatIndx ~= 0;
    xCoord = featuresInfo(iFrame).xCoord(:,1);
    yCoord = featuresInfo(iFrame).yCoord(:,1);
    amp = featuresInfo(iFrame).amp(:,1);
    sX = featuresInfo(iFrame).stdAlong(:,1);
    sY = featuresInfo(iFrame).stdAside(:,1);
    theta = featuresInfo(iFrame).theta(:,1);
    bkg = featuresInfo(iFrame).bkg(:,1);
    
    allTrackParams(ind,1) = xCoord(tracksFeatIndx(ind));
    allTrackParams(ind,2) = yCoord(tracksFeatIndx(ind));
    allTrackParams(ind,3) = amp(tracksFeatIndx(ind));
    allTrackParams(ind,4) = sX(tracksFeatIndx(ind));
    allTrackParams(ind,5) = sY(tracksFeatIndx(ind));
    allTrackParams(ind,6) = theta(tracksFeatIndx(ind));
    allTrackParams(ind,7) = bkg(tracksFeatIndx(ind));
end

%% Interpolate feature parameters in gaps
gacombIdx = diff(isnan(allTrackParams(:,1)));
gapStarts = find(gacombIdx==1)+1;
gapEnds = find(gacombIdx==-1);
gapLengths = gapEnds-gapStarts+1;
nGaps = length(gapLengths);

for g = 1:nGaps
    borderIdx = [gapStarts(g)-1 gapEnds(g)+1];
    gacombIdx = gapStarts(g):gapEnds(g);
    allTrackParams(gacombIdx,1) = interp1(borderIdx, allTrackParams(borderIdx,1), gacombIdx);
    allTrackParams(gacombIdx,2) = interp1(borderIdx, allTrackParams(borderIdx,2), gacombIdx);
    allTrackParams(gacombIdx,3) = interp1(borderIdx, allTrackParams(borderIdx,3), gacombIdx);
    allTrackParams(gacombIdx,4) = interp1(borderIdx, allTrackParams(borderIdx,4), gacombIdx);
    allTrackParams(gacombIdx,5) = interp1(borderIdx, allTrackParams(borderIdx,5), gacombIdx);
    allTrackParams(gacombIdx,6) = interp1(borderIdx, allTrackParams(borderIdx,6), gacombIdx);
    allTrackParams(gacombIdx,7) = interp1(borderIdx, allTrackParams(borderIdx,7), gacombIdx);
end

%% Find the set of track pairs that overlap in time

E = pcombs(1:nTracks);

tOverlapFirst = max(tFirst(E(:,1)), tFirst(E(:,2)));
tOverlapLast = min(tLast(E(:,1)), tLast(E(:,2)));
overlap = tOverlapLast - tOverlapFirst + 1 + timeGap;

hasOverlap = overlap >= minOverlap;

% trim arrays
E = E(hasOverlap,:);
tOverlapFirst = tOverlapFirst(hasOverlap);
tOverlapLast = tOverlapLast(hasOverlap);
overlap = overlap(hasOverlap);

%% Compute the euclidian distance between pair of tracks

% pFirst and pLast are indexing every variable named 'allTrack*'
pLast = cumsum(lifetime);
pFirst = pLast-lifetime+1;

% Point the location of each track parameters for every pair
pFirst1 = pFirst(E(:,1)) + tOverlapFirst - tFirst(E(:,1));
pFirst2 = pFirst(E(:,2)) + tOverlapFirst - tFirst(E(:,2));

% allPairParams1 = arrayfun(@(a,b) allTrackParams(a:a+b-1,1:2),...
%     pFirst1 + tOverlapFirst - tFirst(E(:,1)), overlap,...
%     'UniformOutput', false);
% allPairParams1 = vertcat(allPairParams1{:});
% allPairParams1 = num2cell(allPairParams1,1);
% 
% allPairParams2 = arrayfun(@(a,b) allTrackParams(a:a+b-1,1:2),...
%     pFirst2 + tOverlapFirst - tFirst(E(:,2)), overlap,...
%     'UniformOutput', false);
% allPairParams2 = vertcat(allPairParams2{:});
% allPairParams2 = num2cell(allPairParams2,1);
% 
% % ppFirst and ppLast are indexing every variable named 'allPair*'
% ppLast = cumsum(overlap);
% ppFirst = ppLast-overlap+1;
% 
% [x1,y1] = allPairParams1{:};
% [x2,y2] = allPairParams2{:};
% 
% allPWD = sqrt((x1-x2).^2 + (y1-y2).^2);
% 
% PWD = arrayfun(@(a,b) mean(allPWD(a:a+b-1)), ppFirst, overlap);

% sort overlap values
[overlap idx] = sort(overlap);

% sort ANY other vector that contains pair of tracks information
E = E(idx,:);
tOverlapFirst = tOverlapFirst(idx);
tOverlapLast = tOverlapLast(idx);
pFirst1 = pFirst1(idx);
pFirst2 = pFirst2(idx);

firstIdx = find([1; diff(overlap)]);
lastIdx = find([-diff(-overlap); 1]);

overlapValues = unique(overlap);

PWD = zeros(size(overlap));

for k = 1:length(firstIdx)
    % indexes corresponding to overlap value
    range = firstIdx(k):lastIdx(k);
    
    M = repmat((0:overlapValues(k)-1), [length(range) 1]);
    idx1 = repmat(pFirst1(range), [1 overlapValues(k)]) + M;
    idx2 = repmat(pFirst2(range), [1 overlapValues(k)]) + M;
    
    x1 = allTrackParams(idx1,1);
    y1 = allTrackParams(idx1,2);
    x2 = allTrackParams(idx2,1);
    y2 = allTrackParams(idx2,2);
    
    % average distance
    PWD(range) = mean(reshape(sqrt((x1-x2).^2 + (y1-y2).^2), ...
        [length(range) overlapValues(k)]), 2);
end

% Trim arrays
isCloseEnough = PWD <= maxEuclidianDist;
E = E(isCloseEnough,:);
overlap = overlap(isCloseEnough);
tOverlapFirst = tOverlapFirst(isCloseEnough);
tOverlapLast = tOverlapLast(isCloseEnough);
pFirst1 = pFirst1(isCloseEnough);
pFirst2 = pFirst2(isCloseEnough);

%% Remove any pair of tracks that is:
% - within the first 'bandWidth' nanometers of the cell edge
% AND
% - parallel to the cell edge

distToEdgePath = movieData.distanceTransform.directory;
distToEdgeFiles = dir([distToEdgePath filesep '*.mat']);

isTrackInBand = true(nTracks,1);
bandWidth = bandWidth / movieData.pixelSize_nm;

for iFrame = 1:nFrames
    load(fullfile(distToEdgePath, distToEdgeFiles(iFrame).name));
    
    inputIdx = allTrackFrameIdx == iFrame;
    outputIdx = iFrame >= tFirst & iFrame <= tLast;
    
    xi = round(allTrackParams(inputIdx,1));
    yi = round(allTrackParams(inputIdx,2));
    ind = sub2ind(size(distToEdge),yi,xi);
    
    isTrackInBand(outputIdx) = isTrackInBand(outputIdx) & ...
        distToEdge(ind) < bandWidth;
end

isTrackPairInBand = isTrackInBand(E(:,1)) & ...
    isTrackInBand(E(:,2));

nPixels = zeros(size(overlap));
nPixelsAlongFlow = zeros(size(overlap));

% angular accuracy
p = .5;
dt = p * pi / 2;

for iFrame = 1:nFrames
    % Load distance transform
    load(fullfile(distToEdgePath, distToEdgeFiles(iFrame).name));

    isPairInFrame = iFrame >= tOverlapFirst & iFrame <= tOverlapFirst + overlap - 1;
    isValidPair = isTrackPairInBand & isPairInFrame;

    % Get the coordinates of the extremities of the valid pair
    offset = iFrame - tOverlapFirst(isValidPair);
    ind1 = pFirst1(isValidPair) + offset;
    ind2 = pFirst2(isValidPair) + offset;
    
    x1 = allTrackParams(ind1,1);
    x2 = allTrackParams(ind2,1);
    y1 = allTrackParams(ind1,2);
    y2 = allTrackParams(ind2,2);
    
    % Compute the unit vector along the pair
    dx = x2 - x1;
    dy = y2 - y1;
    norm = sqrt(dx.^2 + dy.^2);
    dx = dx ./ norm;
    dy = dy ./ norm;
    
    % Integer-valued extremity points
    x1i = round(x1);
    y1i = round(y1);
    x2i = round(x2);
    y2i = round(y2);
    
    % Compute Bresenham line between the 2 extremity points
    ptsLines = arrayfun(@(xx1,yy1,xx2,yy2) bresenhamMEX([xx1,yy1], [xx2,yy2]), ...
        x1i, y1i, x2i, y2i, 'UniformOutput', false);

    nPts = cellfun(@(p) size(p,1), ptsLines);
    ptsLast = cumsum(nPts);
    ptsFirst = ptsLast - nPts + 1;
    ptsLines = vertcat(ptsLines{:});

    indLines = sub2ind(imSize, ptsLines(:,2), ptsLines(:,1));

    % Get cell edge orientation along the Bresenham lines
    [dU dV] = gradient(distToEdge);
    dU = dU(indLines);
    dV = dV(indLines);
    norm = sqrt(dU.^2 + dV.^2);
    isnnz = norm ~= 0;
    dU = dU(isnnz) ./ norm(isnnz);
    dV = dV(isnnz) ./ norm(isnnz);
    
    nValidPts = arrayfun(@(a,b) nnz(isnnz(a:a+b-1)), ptsFirst, nPts);
    ptsLast = cumsum(nValidPts);
    ptsFirst = ptsLast - nValidPts + 1;
    
    % Expand dx and dy
    ind = arrayfun(@(a,b) ones(a,1) * b, nValidPts, (1:numel(dx))', ...
        'UniformOutput', false);
    ind = vertcat(ind{:});
    dx = dx(ind);
    dy = dy(ind);
        
    thetaLines = acos(abs(dx .* dU + dy .* dV));
    
    isAlongFlow = thetaLines <= dt;
    
    nPixels(isValidPair) = nPixels(isValidPair) + nValidPts;
    nPixelsAlongFlow(isValidPair) = nPixelsAlongFlow(isValidPair) + ...
        arrayfun(@(a,b) nnz(isAlongFlow(a:a+b-1)), ptsFirst, nValidPts);
end

maxNPixels = max(nPixels);

cutoffs = icdf('bino', 1-alpha, (1:maxNPixels)', p);

isValidPair = true(size(overlap));
isValidPair(isTrackPairInBand) = nPixelsAlongFlow(isTrackPairInBand) >= ...
    cutoffs(nPixels(isTrackPairInBand));

% trim arrays
E = E(isValidPair,:);
tOverlapFirst = tOverlapFirst(isValidPair);
tOverlapLast = tOverlapLast(isValidPair);
overlap = overlap(isValidPair);
pFirst1 = pFirst1(isValidPair);
pFirst2 = pFirst2(isValidPair);

%% DEBUG
segments = cell(nFrames,1);
for iFrame = 1:nFrames
    isPairInFrame = iFrame >= tOverlapFirst & iFrame <= tOverlapFirst + overlap - 1;
    
    % Get the coordinates of the extremities of the valid pair
    offset = iFrame - tOverlapFirst(isPairInFrame);
    ind1 = pFirst1(isPairInFrame) + offset;
    ind2 = pFirst2(isPairInFrame) + offset;
    
    x1 = allTrackParams(ind1,1);
    x2 = allTrackParams(ind2,1);
    y1 = allTrackParams(ind1,2);
    y2 = allTrackParams(ind2,2);
    
    segments{iFrame} = [x1 y1 x2 y2];
end

save(fullfile(movieData.pairTracks.directory, 'pairTrackCands.mat'), 'segments');

%% Compute the edge weight

W = zeros(size(overlap));

allPairParams1 = arrayfun(@(a,b) allTrackParams(a:a+b-1,:),...
    pFirst1 + tOverlapFirst - tFirst(E(:,1)), overlap,...
    'UniformOutput', false);
allPairParams1 = vertcat(allPairParams1{:});
allPairParams1 = num2cell(allPairParams1,1);

allPairParams2 = arrayfun(@(a,b) allTrackParams(a:a+b-1,:),...
    pFirst2 + tOverlapFirst - tFirst(E(:,2)), overlap,...
    'UniformOutput', false);
allPairParams2 = vertcat(allPairParams2{:});
allPairParams2 = num2cell(allPairParams2,1);
[x1,y1,~,sx1,~,t1,~] = allPairParams1{:};
[x2,y2,~,sx2,~,t2,~] = allPairParams2{:};

% WD is the weight associated with pairwise distance between pair of tracks
allPWD = sqrt((x1-x2).^2 + (y1-y2).^2);
WD = exp(-allPWD);

% WA is the weight associated with the deviation of each feature from
% the pair axis
ct = cos(t1);
st = sin(t1);
x0 = x1 - sx1 .* ct;
y0 = y1 - sx1 .* st;
pD = abs(((x2-x1) .* (y1 - y0) - (x1 - x0) .* (y2 - y1)) ./ allPWD);
cos1 = sqrt(1 - (pD ./ sx1).^2);

ct = cos(t2);
st = sin(t2);
x0 = x2 - sx2 .* ct;
y0 = y2 - sx2 .* st;
pD = abs(((x2-x1) .* (y1 - y0) - (x1 - x0) .* (y2 - y1)) ./ allPWD);
cos2 = sqrt(1 - (pD ./ sx2).^2);

WA = .5 * (cos1 + cos2);

ppLast = cumsum(overlap);
ppFirst = ppLast-overlap+1;

W = arrayfun(@(a,b) sum(WD(a:a+b-1) .* WA(a:a+b-1)), ppFirst, overlap);

%% Compute the edge weight
% 
% W = zeros(size(overlap));
% 
% % Compute the segment parameters for each pair of tracks
% allPairParams1 = arrayfun(@(a,b) allTrackParams(a:a+b-1,:),...
%     pFirst1 + tOverlapFirst - tFirst(E(:,1)), overlap,...
%     'UniformOutput', false);
% allPairParams1 = vertcat(allPairParams1{:});
% allPairParams1 = num2cell(allPairParams1,1);
% 
% allPairParams2 = arrayfun(@(a,b) allTrackParams(a:a+b-1,:),...
%     pFirst2 + tOverlapFirst - tFirst(E(:,2)), overlap,...
%     'UniformOutput', false);
% allPairParams2 = vertcat(allPairParams2{:});
% allPairParams2 = num2cell(allPairParams2,1);
% 
% [x1,y1,amp1,sx1,sy1,~,c1] = allPairParams1{:};
% [x2,y2,amp2,sx2,sy2,~,c2] = allPairParams2{:};
% 
% % Compute the support of each track pair
% tC = atan2(y2-y1,x2-x1);
% ct = cos(tC);
% st = sin(tC);
% ex1 = x1 - sx1 .* ct;
% ey1 = y1 - sx1 .* st;
% ex2 = x2 + sx2 .* ct;
% ey2 = y2 + sx2 .* st;
% xC = .5 * (ex1+ex2);
% yC = .5 * (ey1+ey2);
% lC = sqrt((ex2-ex1).^2 + (ey2-ey1).^2);
% sC = .5 * (sy1+sy2);
% cC = .5 * (c1 + c2);
% amp1 = amp1 + (c1 - cC);
% amp2 = amp2 + (c2 - cC);
% 
% [xRangeAll, yRangeAll, nzIdxAll] = arrayfun(@(x0,y0,l0,s0,t0) ...
%     segment2DSupport(x0,y0,l0,s0,t0,kSigma,[nx ny]),...
%     xC,yC,lC,sC,tC,'UniformOutput', false);
% 
% % Compute the error variance of each pair model
% imagePath = fullfile(movieData.imageDirectory, ...
%     movieData.channelDirectory{iChannel});
% imageFiles = dir([imagePath filesep '*.tif']);
% 
% for iFrame = 1:nFrames
%     % Read image
%     ima = double(imread(fullfile(imagePath, imageFiles(iFrame).name)));
%         
%     % Find which pair of CC is living at frame iFrame
%     isPairInFrame = iFrame >= tOverlapFirst & iFrame <= tOverlapFirst + overlap - 1;
%     
%     % Compute the indices of arrays allPairParams*, xRange, yRange, nzIdx
%     ppIdx = ppFirst(isPairInFrame) + iFrame - tOverlapFirst(isPairInFrame);
%     
%     % Crop the image for the first feature of each pair
%     xRange = xRangeAll(ppIdx);
%     yRange = yRangeAll(ppIdx);
%     nzIdx = nzIdxAll(ppIdx);
%     
%     imaCrops = cellfun(@(xR,yR) ima(yR,xR), xRange, yRange,...
%         'UniformOutput', false);
%     allPixels = cellfun(@(I,p) I(p), imaCrops, nzIdx,...
%         'UniformOutput', false);
% 
%     % Concatenate every model 2 values (pair model)
%     allModel1 = arrayfun(@(x,y,amp,sx,sy,t,xR,yR,nzIdx)...
%         anisoGaussian2D(x,y,amp,sx,sy,t,xR{1},yR{1},nzIdx{1}),...
%         x1(ppIdx), y1(ppIdx), amp1(ppIdx), sx1(ppIdx), sy1(ppIdx),...
%         tC(ppIdx), xRange, yRange, nzIdx,'UniformOutput', false);
%     
%     allModel2 = arrayfun(@(x,y,amp,sx,sy,t,xR,yR,nzIdx)...
%         anisoGaussian2D(x,y,amp,sx,sy,t,xR{1},yR{1},nzIdx{1}),...
%         x2(ppIdx), y2(ppIdx), amp2(ppIdx), sx2(ppIdx), sy2(ppIdx),...
%         tC(ppIdx), xRange, yRange, nzIdx,'UniformOutput', false);
%     
%     vInv = arrayfun(@(I,m1,m2,bkg) ...
%         (numel(I{1}) - 1) / (sum((I{1} - m1{1} - m2{1} - bkg).^2)),...
%         allPixels, allModel1, allModel2, cC(ppIdx));
%     
%     M = maxWeightedMatching(nTracks,E(isPairInFrame,:), vInv);
%     
%     W(isPairInFrame) = W(isPairInFrame) + vInv .* M;
% end

%% Pairing

% Compute Maximum Weighted Matching
M = maxWeightedMatching(nTracks,E,W);

% Trim
E = E(M,:);
tOverlapFirst = tOverlapFirst(M);
tOverlapLast = tOverlapLast(M);
overlap = overlap(M);
pFirst1 = pFirst1(M);
pFirst2 = pFirst2(M);

%% DEBUG
segments = cell(nFrames,1);
for iFrame = 1:nFrames
    isPairInFrame = iFrame >= tOverlapFirst & iFrame <= tOverlapFirst + overlap - 1;
    
    % Get the coordinates of the extremities of the valid pair
    offset = iFrame - tOverlapFirst(isPairInFrame);
    ind1 = pFirst1(isPairInFrame) + offset;
    ind2 = pFirst2(isPairInFrame) + offset;
    
    x1 = allTrackParams(ind1,1);
    x2 = allTrackParams(ind2,1);
    y1 = allTrackParams(ind1,2);
    y2 = allTrackParams(ind2,2);
    
    segments{iFrame} = [x1 y1 x2 y2];
end

save(fullfile(movieData.pairTracks.directory, 'pairTracks.mat'), 'segments');

%% END
movieData.pairTracks.status = 1;
