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

% pFirst and pLast are indexing every variable named 'all*')
pLast = cumsum(lifetime);
pFirst = pLast-lifetime+1;

% Point the location of each track parameters for every pair
pFirst1 = pFirst(E(:,1)) + tOverlapFirst - tFirst(E(:,1));
pFirst2 = pFirst(E(:,2)) + tOverlapFirst - tFirst(E(:,2));

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

euclidianDist = zeros(size(overlap));

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
    euclidianDist(range) = mean(reshape(sqrt((x1-x2).^2 + (y1-y2).^2), ...
        [length(range) overlapValues(k)]), 2);
end

% Trim arrays
isCloseEnough = euclidianDist <= maxEuclidianDist;
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

% %% DEBUG
% segments = cell(nFrames,1);
% for iFrame = 1:nFrames
%     isPairInFrame = iFrame >= tOverlapFirst & iFrame <= tOverlapFirst + overlap - 1;
%     
%     % Get the coordinates of the extremities of the valid pair
%     offset = iFrame - tOverlapFirst(isPairInFrame);
%     ind1 = pFirst1(isPairInFrame) + offset;
%     ind2 = pFirst2(isPairInFrame) + offset;
%     
%     x1 = allTrackParams(ind1,1);
%     x2 = allTrackParams(ind2,1);
%     y1 = allTrackParams(ind1,2);
%     y2 = allTrackParams(ind2,2);
%     
%     segments{iFrame} = [x1 y1 x2 y2];
% end
% 
% save(fullfile(movieData.pairTracks.directory, 'CCpairCands.mat'), 'segments');

%% Compute the residue of model 1 and 2 for each pair of tracks

% Compute the support of each track feature
cellParams = num2cell(allTrackParams(:,[1 2 4 5 6]),1);
[xRange, yRange, nzIdx] = arrayfun(@(x0,y0,sX,sY,theta) ...
    modelSupportFunc{1}(x0,y0,sX,sY,theta,kSigma,[nx ny]),...
    cellParams{:},'UniformOutput', false);

% Compute the segment parameters for each pair of tracks
allSegmentParams1 = arrayfun(@(a,b) allTrackParams(a:a+b-1,:),...
    pFirst1 + tOverlapFirst - tFirst(E(:,1)), overlap,...
    'UniformOutput', false);
allSegmentParams1 = vertcat(allSegmentParams1{:});
allSegmentParams1 = num2cell(allSegmentParams1,1);

allSegmentParams2 = arrayfun(@(a,b) allTrackParams(a:a+b-1,:),...
    pFirst2 + tOverlapFirst - tFirst(E(:,2)), overlap,...
    'UniformOutput', false);
allSegmentParams2 = vertcat(allSegmentParams2{:});
allSegmentParams2 = num2cell(allSegmentParams2,1);

[x1,y1,amp1,sx1,sy1,t1,bkg1] = allSegmentParams1{:};
[x2,y2,amp2,sx2,sy2,t2,bkg2] = allSegmentParams2{:};

tC = atan2(y2-y1,x2-x1);

% Compute the residue of each pair against model 1 & 2
imagePath = fullfile(movieData.imageDirectory, ...
    movieData.channelDirectory{iChannel});
imageFiles = dir([imagePath filesep '*.tif']);

% ppFirst is an array nPairs X 1 pointing to the first element of each pair
% in allSegmentParams*
ppLast = cumsum(overlap);
ppFirst = ppLast-overlap+1;

% residue for model 1 and 2 for each pair
R1 = zeros(size(overlap));
R2 = zeros(size(overlap));
% number of pixels per pair (segment support over the overlap)
N = zeros(size(overlap));

for iFrame = 1:nFrames
    % Read image
    ima = double(imread(fullfile(imagePath, imageFiles(iFrame).name)));
        
    % Find which pair of CC is living at frame iFrame
    isPairInFrame = iFrame >= tOverlapFirst & iFrame <= tOverlapFirst + overlap - 1;
    
    % Compute the indeices of arrays xRange, yRange, nzIdx
    pIdx1 = pFirst1(isPairInFrame) + iFrame - tOverlapFirst(isPairInFrame);
    pIdx2 = pFirst2(isPairInFrame) + iFrame - tOverlapFirst(isPairInFrame);

    % Compute the indices of arrays allSegmentParams*
    ppIdx = ppFirst(isPairInFrame) + iFrame - tOverlapFirst(isPairInFrame);
    
    % Crop the image for the first feature of each pair
    xRange1 = xRange(pIdx1);
    yRange1 = yRange(pIdx1);
    nzIdx1 = nzIdx(pIdx1);

    imaCrops = cellfun(@(xR,yR) ima(yR,xR), xRange1, yRange1,...
        'UniformOutput', false);
    allPixels1 = cellfun(@(I,p) I(p), imaCrops, nzIdx1,...
        'UniformOutput', false);

    % Crop the image for the 2nd featue of each pair
    xRange2 = xRange(pIdx2);
    yRange2 = yRange(pIdx2);
    nzIdx2 = nzIdx(pIdx2);
    
    imaCrops = cellfun(@(xR,yR) ima(yR,xR), xRange2, yRange2,...
        'UniformOutput', false);
    allPixels2 = cellfun(@(I,p) I(p), imaCrops, nzIdx2,...
        'UniformOutput', false);
    
    % Compute model 1 (individual tracks model)
    allModel11 = arrayfun(@(x,y,amp,sx,sy,t,bkg,xR,yR,nzIdx)...
        anisoGaussian2D(x,y,amp,sx,sy,t,xR{1},yR{1},nzIdx{1}) + bkg,...
        x1(ppIdx), y1(ppIdx), amp1(ppIdx), sx1(ppIdx), sy1(ppIdx),...
        t1(ppIdx), bkg1(ppIdx), xRange1, yRange1, nzIdx1,...
        'UniformOutput', false);
    
    allModel12 = arrayfun(@(x,y,amp,sx,sy,t,bkg,xR,yR,nzIdx)...
        anisoGaussian2D(x,y,amp,sx,sy,t,xR{1},yR{1},nzIdx{1}),...
        x2(ppIdx), y2(ppIdx), amp2(ppIdx), sx2(ppIdx), sy2(ppIdx),...
        t2(ppIdx), bkg2(ppIdx), xRange2, yRange2, nzIdx2,...
        'UniformOutput', false);
    
    % Concatenate every model 2 values (pair model)
    allModel21 = arrayfun(@(x,y,amp,sx,sy,t,bkg,xR,yR,nzIdx)...
        anisoGaussian2D(x,y,amp,sx,sy,t,xR{1},yR{1},nzIdx{1}) + bkg,...
        x1(ppIdx), y1(ppIdx), amp1(ppIdx), sx1(ppIdx), sy1(ppIdx),...
        tC(ppIdx), bkg1(ppIdx), xRange1, yRange1, nzIdx1,...
        'UniformOutput', false);
    
    allModel22 = arrayfun(@(x,y,amp,sx,sy,t,bkg,xR,yR,nzIdx)...
        anisoGaussian2D(x,y,amp,sx,sy,t,xR{1},yR{1},nzIdx{1}),...
        x2(ppIdx), y2(ppIdx), amp2(ppIdx), sx2(ppIdx), sy2(ppIdx),...
        tC(ppIdx), bkg2(ppIdx), xRange2, yRange2, nzIdx2,...
        'UniformOutput', false);
    
    % Compute the residue of model 1 and model 2
    R1(isPairInFrame) = R1(isPairInFrame) + cellfun(@(I1,I2,m11,m12) ...
        sum((I1 - m11).^2) + sum((I2 - m12).^2), ...
        allPixels1, allPixels2, allModel11, allModel12);
    
    R2(isPairInFrame) = R2(isPairInFrame) + cellfun(@(I1,I2,m21,m22) ...
        sum((I1 - m21).^2) + sum((I2 - m22).^2), ...
        allPixels1, allPixels2, allModel21, allModel22);
    
    N(isPairInFrame) = N(isPairInFrame) + cellfun(@numel, nzIdx1) + ...
        cellfun(@numel, nzIdx2);
end

R1 = R1 ./ (N - 1);
R2 = R2 ./ (N - 1);

%% Remove pair where model 2 is significantly better than model 1
assert(all(R1 < R2));


% % CC: a cell array of connected components. Each connected component
% % contains a list of track indices. At the beginning, Each connected
% % component contains one track only.
% CC = arrayfun(@(x) {x}, (1:nTracks)');
% nCC = numel(CC);
% 
% % CCmodels: is a struct array containing the following fields:
% %
% % .tFirst:              array of nCCx1 corresponding to the first common
% %                       frame of every track belonging to the CC
% %
% % .tLast:               array of nCCx1 corresponding to the last common
% %                       frame of every track belonging to the CC
% %
% % .type:                array of nCCx1 taking the following values:
% %                       1 == single anisotropic feature
% %                       2 == pair of anisotropic features
% %
% % .pFirst:              array of nCCx1 containing the first index in
% %                       .allModelParams{*} array. The index in the cell
% %                       array .allModelParams{*} is given by .modelType.
% %
% % .allModelParams:      cell array of size 2 containing model parameters of
% %                       anisotropic models (index 1) and segment models
% %                       (index 2):
% %
% %                       if modelType == 1:
% %                          x1 y1
% %                          x2 y2
% %                          ...
% %                          xN yN
% %
% %                       if modelType == 2:
% %                          x1 y1 l1 T1
% %                          x2 y2 l2 T2
% %                          ...
% %                          xN yN lN TN
% %
% % .xRange:              cell array of size 2 containing model support
% %                       (xRange) of anisotropic models (index 1) and
% %                       segment models (index 2)
% %
% % .yRange:              cell array of size 2 containing model support
% %                       (yRange) of anisotropic models (index 1) and
% %                       segment models (index 2)
% %
% % .nzIdx:               cell array of size 2 containing model support
% %                       (nzIdx) of anisotropic models (index 1) and segment
% %                       models (index 2)
% 
% CCmodels.tFirst = tFirst;
% CCmodels.tLast = tLast;
% CCmodels.type = ones(1, nCC);
% CCmodels.pFirst = first;
% CCmodels.allModelParams = cell(2,1);
% CCmodels.allModelParams{1} = allTrackParams(1:2,:);
% 
% % Compute the model support of each CC
% CCmodels.xRange = cell(2,1);
% CCmodels.yRange = cell(2,1);
% CCmodels.nzIdx = cell(2,1);
% 
% cellParams = num2cell(allTrackParams(:,[1 2 4 5 6]),1);
% [xRange, yRange, nzIdx] = arrayfun(@(x0,y0,sX,sY,theta) ...
%     modelSupportFunc{1}(x0,y0,sX,sY,theta,kSigma,[nx ny]),...
%     cellParams{:},'UniformOutput', false);
% CCmodels.xRange{1} = xRange;
% CCmodels.yRange{1} = yRange;
% CCmodels.nzIdx{1} = nzIdx;
% 
% % pairCCIdx: pair of CC candidates. We initialize pairCCIdx with E
% pairCCIdx = E;
% nPairCC = size(pairCCIdx,1);
% 
% for iLevel = 1:nLevels
% 
%     % Compute model parameters of each pair of CC
%     iFirst1 = CCmodels.tFirst(pairCCIdx(:,1));
%     iFirst2 = CCmodels.tFirst(pairCCIdx(:,2));
%     iLast1 = CCmodels.tLast(pairCCIdx(:,1));
%     iLast2 = CCmodels.tLast(pairCCIdx(:,2));
%     
%     % pairCCmodels has EXACTLY the same structure as CCmodels but instead
%     % of representing models of CC, it represents models of pairs of CC.
%     pairCCmodels.tFirst = max(iFirst1,iFirst2);
%     pairCCmodels.tLast = min(iLast1,iLast2);
%     pairCCmodels.type = repmat(2,1,nPairCC);    
%     overlap = pairCCmodels.tLast - pairCCmodels.tFirst + 1;
%     offset1 = pairCCmodels.tFirst - iFirst1;
%     offset2 = pairCCmodels.tFirst - iFirst2;        
%     pairCCmodels.pFirst = cumsum([1 overlap(1:end-1)]);
%     pairCCmodels.params = cell(2,1);
%     pairCCmodels.params{2} = zeros(sum(overlap),7);
%     
%     % iFunc == 1: modelType1 == 1, modelType2 == 1
%     % iFunc == 2: modelType1 == 1, modelType2 == 2
%     % iFunc == 3: modelType1 == 2, modelType2 == 1
%     % iFunc == 4: modelType1 == 2, modelType2 == 2
%     for iFunc = 1:4
%         modelType1 = ceil(iFunc/2);
%         modelType2 = mod(iFunc-1,2)+1;
%         
%         isFunc = ...
%             CCmodels.type(pairCCIdx(:,1)) == modelType1 & ...
%             CCmodels.type(pairCCIdx(:,2)) == modelType2;
%         
%         inputInd1 = cell2mat(arrayfun(@(x,y) x:(x+y-1), ...
%             CCmodels.pFirst(pairCCIdx(isFunc,1)) + offset1(isFunc), ...
%             overlap(isFunc), 'UniformOutput',false));
%         
%         inputInd2 = cell2mat(arrayfun(@(x,y) x:(x+y-1), ...
%             CCmodels.pFirst(pairCCIdx(isFunc,2)) + offset2(isFunc), ...
%             overlap(isFunc), 'UniformOutput', false));
%         
%         outputInd = cell2mat(arrayfun(@(x,y) x:(x+y-1), ...
%             pairCCmodels.pFirst(isFunc), overlap(isFunc), ...
%             'UniformOutput', false));
%         
%         pairCCmodels.params{2}(outputInd,:) = initFunc{iFunc}(...
%             CCmodels.params{modelType1}(inputInd1,:), ...
%             CCmodels.params{modelType2}(inputInd2,:));
%     end
%     
% %     if ~batchMode
% %         h = waitbar(0,['Please wait, computing pair models at level' ...
% %             num2str(iLevel) '...']);
% %     end
% % 
% %     W = zeros(1,nPairCC);
% %     
% %     % Compute pair model
% %     for iFrame = 1:nFrames
% %         % Read image
% %         ima = double(imread(fullfile(imagePath, imageFiles(iFrame).name)));
% %         
% %         % Find which pair of CC is living at frame iFrame
% %         pairInFrame = find(pairCCmodels.tFirst <= iFrame & pairCCmodels.tLast >= iFrame);
% % 
% %         % some variables on individual CCs associated by the pair candidate
% %         cc1 = pairCCIdx(pairInFrame,1);
% %         cc2 = pairCCIdx(pairInFrame,2);
% %         p1Idx = CCmodels.pFirst(cc1) + iFrame - CCmodels.tFirst(cc1);
% %         p2Idx = CCmodels.pFirst(cc2) + iFrame - CCmodels.tFirst(cc2);
% %         modelType1 = CCmodels.type(cc1);
% %         modelType2 = CCmodels.type(cc2);
% %         
% %         pIdx = pairCCmodels.pFirst(pairInFrame) + iFrame - pairCCmodels.tFirst(pairInFrame);
%         
% %         initParams = pairCCmodels.params{2}(pIdx,:);
% %         nInitParams = size(initParams);
%         
% %         for iP = 1:nInitParams
% %             
% %             % Fit model
% %             cellParams = num2cell(initParams(iP,:),1);
% %             [x0 y0 A l sigma theta C] = cellParams{:};
%             
% %             [xRange,yRange,nzIdx] = modelSupportFunc{2}(x0,y0,l,sigma,theta,kSigma,[nx ny]);
% 
% %             crop = ima(yRange,xRange);
% %             mask = false(size(crop));
% %             mask(nzIdx) = true;
% %             crop(~mask) = NaN;            
%             
% %             xi = xRange(floor(numel(xRange)/2));
% %             yi = yRange(floor(numel(yRange)/2));            
% %             prm = fitModelFunc{2}(crop, [x0-xi,y0-yi,A,l,sigma,theta,C], 'AlC');
% %             
% %             % Assign new fitted parameters
% %             prm(1:2) = [xi, yi] + prm(1:2);            
% %             cellParams = num2cell(prm,1);
% %             pairCCmodels.params{2}(pIdx(iP),:) = prm;
%             
% %             % Get the support of CC 1
% %             xRange1 = CCmodels.xRange{modelType1(iP)}{p1Idx(iP)};
% %             yRange1 = CCmodels.yRange{modelType1(iP)}{p1Idx(iP)};
% %             nzIdx1 = CCmodels.nzIdx{modelType1(iP)}{p1Idx(iP)};
% %             mask1 = false(numel(yRange1),numel(xRange1));
% %             mask1(nzIdx1) = true;
% %             
% %             % Get the support of CC 2
% %             xRange2 = CCmodels.xRange{modelType2(iP)}{p2Idx(iP)};
% %             yRange2 = CCmodels.yRange{modelType2(iP)}{p2Idx(iP)};
% %             nzIdx2 = CCmodels.nzIdx{modelType2(iP)}{p2Idx(iP)};
% %             mask2 = false(numel(yRange2),numel(xRange2));
% %             mask2(nzIdx2) = true;
% %             
% %             % Get the union of CC1 and CC2 support and the support of the
% %             % pair between CC1 and CC2.
% %             xRangeFull = min([xRange(1),xRange1(1),xRange2(1)]):...
% %                 max([xRange(end),xRange1(end),xRange2(end)]);
% %             yRangeFull = min([yRange(1),yRange1(1),yRange2(1)]):...
% %                 max([yRange(end),yRange1(end),yRange2(end)]);
% %         
% %             maskFull = false(numel(yRangeFull),numel(xRangeFull));
% %             tmp1 = yRange1 - yRangeFull(1) + 1;
% %             tmp2 = xRange1 - xRangeFull(1) + 1;
% %             maskFull(tmp1,tmp2) = maskFull(tmp1,tmp2) | mask1;
% %             
% %             tmp1 = yRange2 - yRangeFull(1) + 1;
% %             tmp2 = xRange2 - xRangeFull(1) + 1;
% %             maskFull(tmp1,tmp2) = maskFull(tmp1,tmp2) | mask2;
% %             
% %             tmp1 = yRange - yRangeFull(1) + 1;
% %             tmp2 = xRange - xRangeFull(1) + 1;
% %             maskFull(tmp1,tmp2) = maskFull(tmp1,tmp2) | mask;
% %             
% %             nzIdxFull = find(maskFull);
% %             nPixels = numel(nzIdxFull);
% %             
% %             % Evaluate models of CC1, CC2 and pair between CC1 and CC2
% %             cellParamsCC1 = num2cell(CCmodels.params{modelType1(iP)}(p1Idx(iP),:),1);
% %             cellParamsCC2 = num2cell(CCmodels.params{modelType2(iP)}(p2Idx(iP),:),1);
% %             modelCC1 = modelFunc{modelType1(iP)}(cellParamsCC1{1:end-1},xRangeFull,yRangeFull,nzIdxFull);
% %             modelCC2 = modelFunc{modelType2(iP)}(cellParamsCC2{1:end-1},xRangeFull,yRangeFull,nzIdxFull);
% %             modelPair = modelFunc{2}(cellParams{1:end-1},xRangeFull,yRangeFull,nzIdxFull);
% %             
% %             % Add the background values
% %             modelCC1CC2 = modelCC1 + modelCC2 + ...
% %                 .5 * (cellParamsCC1{end} + cellParamsCC2{end});
% %             
% %             modelPair = modelPair + cellParams{end};
% %             
% %             cropFull = ima(yRangeFull,xRangeFull);
% %             
% %             % Evaluate error variance
% %             varPair = (1/(nPixels-1)) * ...
% %                 sum((cropFull(nzIdxFull) - modelPair(nzIdxFull)).^2);
% % 
% %             varCC1CC2 = (1/(nPixels-1)) * ...
% %                 sum((cropFull(nzIdxFull) - modelCC1CC2(nzIdxFull)).^2);
% %             
% %             % Compute the weight
% %             W(pairInFrame(iP)) = W(pairInFrame(iP)) + ...
% %                 nPixels * log(varCC1CC2 / varPair) + 6 * log(nPixels);
% %         end
%         
% %         if ~batchMode && ishandle(h)
% %             waitbar(iFrame/nFrames, h)
% %         end
% %     end
%     
% %     % Average the weight over the overlap frames
% %     W = bsxfun(@rdivide,W,overlap);
% %     
% %     % Remove pairs with a negative weights
% %     validWeights = W > 0;
%     
% %     if ~batchMode && ishandle(h)
% %         close(h);
% %     end
% 
%     % DEBUG: save pairCCmodels as segment2D
%     segments = cell(nFrames,1);
%     
%     for iFrame = 1:nFrames
%         pairInFrame = find(pairCCmodels.tFirst <= iFrame & ...
%             pairCCmodels.tLast >= iFrame);% & validWeights);
%         
%         pIdx = pairCCmodels.pFirst(pairInFrame) + iFrame - ...
%             pairCCmodels.tFirst(pairInFrame);
%         
%         segments{iFrame} = pairCCmodels.params{2}(pIdx,:);
%     end
%     
%     save(fullfile(movieData.pairTracks.directory, 'CCpairCands.mat'), 'segments');
%     
%     % Compute Maximum Matching
%     D = sparse(pairCCIdx(:,1),pairCCIdx(:,2),W,nCC,nCC,numel(W));
%     M = maxWeightMatching(D,eps);
%     
%     % Accept CCpairParams(M) and merge pairCCIdx(M,1) with pairCCIdx(M,2)
%     % TODO
%     
%     % merge CCparams likewise and replace the merged CC params with
%     % CCPairParams
%     % TODO
%     
%     % Get a new set of CC pairs candidates
%     % TODO
% end

movieData.pairTracks.status = 1;

function pairParams = initPairParams00(params1,params2)
params1 = num2cell(params1,1);
params2 = num2cell(params2,1);
[x1,y1,A1,sx1,sy1,~,C1] = params1{:};
[x2,y2,A2,sx2,sy2,~,C2] = params2{:};

t0 = atan2(y2-y1,x2-x1);
ct = cos(t0);
st = sin(t0);
ex1 = x1 - sx1 .* ct;
ey1 = y1 - sx1 .* st;
ex2 = x2 + sx2 .* ct;
ey2 = y2 + sx2 .* st;
x0 = .5 * (ex1+ex2);
y0 = .5 * (ey1+ey2);
A0 = .5 * (A1+A2);
l0 = sqrt((ex2-ex1).^2 + (ey2-ey1).^2);
s0 = .5 * (sy1+sy2);
C0 = .5 * (C1+C2);
pairParams = [x0 y0 A0 l0 s0 t0 C0];

function pairParams = initPairParams01(params1,params2)
% TODO !!!
pairParams = zeros(size(params1,1),7);

function pairParams = initPairParams11(params1,params2)
% TODO !!!
pairParams = zeros(size(params1,1),7);

% 
% % thetaTracks is a cell array of size nTracks, where each cell contains an
% % array of size = lifetime of the track and contains the local orientaion
% % along the track from each frame
% thetaTracks = arrayfun(@(track) zeros(1, size(track.tracksFeatIndxCG,2)), ...
%     tracksFinal, 'UniformOutput', false);
%     
% for iFrame = 1:nFrames
%     % Get all the tracks that are in iFrame
%     isInFrame = SEL(:,1) <= iFrame & SEL(:,2) >= iFrame;
%     indTrack = find(isInFrame);
%     
%     indFeature = iFrame-SEL(isInFrame,1) + 1;
%     
%     % Get all positions
%     x = X(first(indTrack) + indFeature' - 1);
%     y = Y(first(indTrack) + indFeature' - 1);
%         
%     % Filter image
%     % CHANGE THIS !!!
%     ima = double(imread(fullfile(imagePath, imageFiles(iFrame).name)));
%     [~,T] = steerableFiltering(ima,2,sigmaPSF);
%     ct = cos(T);
%     st = sin(T);
%     
%     % Interpolate vector component at (x,y) positions
%     dU = interp2(ct, x, y,'cubic');
%     dV = interp2(st, x, y,'cubic');
%     
%     % Normalize vector
%     norm = sqrt(dU.^2 + dV.^2);     % norm cannot be 0
%     dU = bsxfun(@rdivide,dU,norm);
%     dV = bsxfun(@rdivide,dV,norm);
%     
%     theta = atan2(dU,dV);
%     
%     % Store theta
%     for iiTrack = 1:numel(indTrack)
%         iTrack = indTrack(iiTrack);
%         thetaTracks{iTrack}(indFeature(iiTrack)) = theta(iiTrack);
%     end
% end
% 
% % compute the cutoff values so that the binomial tail = alpha
% cutoffs = icdf('bino', 1-alpha, 1:nFrames, accT);
% 
% isFeatureAligned = false(nTracks,1);
% 
% wedges = linspace(-pi/2, pi/2 - pi * accT, accT^-1 + 1);
% 
% for iTrack = 1:nTracks
%     theta = thetaTracks{iTrack};
%     
%     isInWedge = bsxfun(@(t,w) t >= w & t < w + pi * accT, theta', wedges);
%     
%     k = max(sum(isInWedge,1));
% 
%     isFeatureAligned(iTrack) = k > cutoffs(numel(theta));
% endD, eps, nonLinkMarker
% 
% % Save tracks classification
% percent = nnz(~isFeatureAligned) * 100 / numel(isFeatureAligned);
% fprintf('Feature-based track classification: static = %f %%, linear = %f %%\n',...
%     percent, 100 - percent);
% 
% trackLabels = isFeatureAligned + 1; %#ok<NASGU> % (1 == static, 2 == aligned)
% save(fullfile(movieData.pairTracks.directory, ...
%     'classifiedTracksPerFeatureAlignement.mat'),'tracksFinal', 'trackLabels');
% 
% % 4) Classify each track whether its dynamics linear or not
% 
% % thetaTracks is a cell array of size nTracks, where each cell contains an
% % array of size = lifetime of the track -1 and contains the direction from
% % one track point to the next
% thetaTracks = cell(nTracks,1);
% 
% for iTrack = 1:nTracks
%     dU = X(first(iTrack)+1:last(iTrack)) - X(first(iTrack):last(iTrack)-1);
%     dV = Y(first(iTrack)+1:last(iTrack)) - Y(first(iTrack):last(iTrack)-1);
%     
%     % Normalize vector
%     norm = sqrt(dU.^2 + dV.^2);     % norm can be 0 !!!
%     
%     theta = zeros(size(norm));
%     isNull = abs(norm) < eps;
%     theta(isNull) = NaN;
%     
%     dU = bsxfun(@rdivide,dU(~isNull),norm(~isNull));
%     dV = bsxfun(@rdivide,dV(~isNull),norm(~isNull));
%     theta(~isNull) = atan2(dU,dV);
%     
%     thetaTracks{iTrack} = theta;
% end
% 
% isDynamicsAligned = false(nTracks,1);
% 
% for iTrack = 1:nTracks
%     theta = thetaTracks{iTrack};
%     
%     nonNaNTheta = theta(~isnan(theta));
% 
%     if ~isempty(nonNaNTheta)
%         isInWedge = bsxfun(@(t,w) t >= w & t < w + pi * accT, nonNaNTheta, wedges);
%     
%         k = max(sum(isInWedge,1));
% 
%         isDynamicsAligned(iTrack) = k > cutoffs(numel(theta));
%     end
% end
% 
% % Save tracks classification
% percent = nnz(~isDynamicsAligned) * 100 / numel(isDynamicsAligned);
% fprintf('Dynamics-based track classification: static = %f %%, linear = %f %%\n',...
%     percent, 100 - percent);
% 
% trackLabels = isDynamicsAligned + 1; %#ok<NASGU> % (1 == static, 2 == aligned)
% save(fullfile(movieData.pairTracks.directory, ...
%     'classifiedTracksPerDynamicsAlignement.mat'),'tracksFinal', 'trackLabels');

