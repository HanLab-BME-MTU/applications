function movieData = getMoviePairTracks(movieData,iChannel,bandWidth,...
    minOverlap,maxEuclidianDist,kSigma,alpha,batchMode)

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
% maxEuclidianDist:     maximum euclidian distance between the mean
%                       position of 2 tracks to be considered as potential
%                       pair of tracks.
%
% kSigma:               cutoff in number of standard deviations
%
% alpha:                quantile of PDF tail. Default is 0.05.

%% Parse input parameters

if nargin < 3 || isempty(bandWidth)
    bandWidth = 1000;
end

if nargin < 4 || isempty(minOverlap)
    minOverlap = 1;
end

if nargin < 7 || isempty(alpha)
    alpha = 0.05;
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

clear kalmanFunctions kalmanInfoLink gapCloseParam costMatrices ...
    seqOfEvents SEL;

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

clear featuresInfo tracksFinal tracksFeatIndx ind xCoord yCoord amp sX ...
    sY theta bkg iFrame;

%% Interpolate feature parameters in gaps
gacombIdx = diff(isnan(allTrackParams(:,1)));
gapStarts = find(gacombIdx==1)+1;
gapEnds = find(gacombIdx==-1);
gapLengths = gapEnds-gapStarts+1;
nGaps = length(gapLengths);

for iGap = 1:nGaps
    borderIdx = [gapStarts(iGap)-1 gapEnds(iGap)+1];
    gacombIdx = gapStarts(iGap):gapEnds(iGap);
    allTrackParams(gacombIdx,1) = interp1(borderIdx, allTrackParams(borderIdx,1), gacombIdx);
    allTrackParams(gacombIdx,2) = interp1(borderIdx, allTrackParams(borderIdx,2), gacombIdx);
    allTrackParams(gacombIdx,3) = interp1(borderIdx, allTrackParams(borderIdx,3), gacombIdx);
    allTrackParams(gacombIdx,4) = interp1(borderIdx, allTrackParams(borderIdx,4), gacombIdx);
    allTrackParams(gacombIdx,5) = interp1(borderIdx, allTrackParams(borderIdx,5), gacombIdx);
    allTrackParams(gacombIdx,6) = interp1(borderIdx, allTrackParams(borderIdx,6), gacombIdx);
    allTrackParams(gacombIdx,7) = interp1(borderIdx, allTrackParams(borderIdx,7), gacombIdx);
end

clear gacombIdx gapStarts gapEnds gapLengths nGaps borderIdx gacombIdx iGap;

%% Find the set of track pairs that overlap in time
E = pcombs(1:nTracks);

tOverlapFirst = max(tFirst(E(:,1)), tFirst(E(:,2)));
tOverlapLast = min(tLast(E(:,1)), tLast(E(:,2)));
overlap = tOverlapLast - tOverlapFirst + 1;
hasOverlap = overlap >= minOverlap;

% trim arrays
E = E(hasOverlap,:);
tOverlapFirst = tOverlapFirst(hasOverlap);
overlap = overlap(hasOverlap);

clear tOverlapLast hasOverlap;

%% Compute the pairwise distance between tracks

% pFirst and pLast are indexing every variable named 'allTrack*'
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
pFirst1 = pFirst1(isCloseEnough);
pFirst2 = pFirst2(isCloseEnough);

clear idx firstIdx lastIdx overlapValues PWD range M idx1 idx2 ...
    x1 y1 x2 y2 isCloseEnough k;

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
overlap = overlap(isValidPair);
pFirst1 = pFirst1(isValidPair);
pFirst2 = pFirst2(isValidPair);

clear distToEdgePath distToEdgeFiles isTrackInBand bandWidth iFrame ...
    allTrackFrameIdx inputIdx outputIdx xi yi ind isTrackPairInBand ...
    nPixels nPixelsAlongFlow p dt isPairInFrame isValidPair offset ...
    ind1 ind2 x1 x2 y1 y2 dx dy norm x1i y1i x2i y2i ptsLines nPts ...
    ptsLast ptsFirst ptsLines indLines dU dV isnnz nValidPts ind...
    thetaLines isAlongFlow maxNPixels cutoffs;

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

clear segments iFrame isPairInFrame offset ind1 ind2 x1 x2 y1 y2;

% IMPORTANT: from here on, we can clear tOverlapFirst, overlap, pFirst1 and
% pFirst2 since they will be recomputed during each next iteration.
clear tOverlapFirst overlap pFirst1 pFirst2;

%% Iterative track clustering

% CC: a cell array of connected components. Each connected component
% contains a list of track indices. At the beginning, Each connected
% component contains one track only.
CC = arrayfun(@(x) {x}, (1:nTracks)');
nCC = numel(CC);

% IMPORTANT: from here on, E points to CC index

for iLevel = 1:1
    
    % The time span of the pair does not cover each track in each CC
    % evenly. Example:
    %
    % CC1 contains tracks 1, 2, 3
    % CC2 contains tracks 4, 5
    %
    % track 1: tFirst = 1, tLast = 5
    % track 2: tFirst = 3, tLast = 9
    % track 3: tFirst = 7, tLast = 14
    % track 4: tFirst = 6, tLast = 20
    % track 5: tFirst = 17, tLast = 30
    %
    % tOverlapFirst = max(min(1,3,7), min(6,17)) = 6
    % tOverlapLast = min(max(5,9,14), max(20,30)) = 14
    % overlap = 9
    %
    % track 1: is not included in the time span of the CC pair
    % track 2: partially included in the time span of the CC pair (6-9)
    % track 3: fully included in the time span of the CC pair
    % track 4: partially included in the time span of the CC pair (6-14)
    % track 5: fully excluded from the time span of the CC pair
    
    % Compute the time span of each CC
    tOverlapFirst = cellfun(@(cc1,cc2) max(min(tFirst(cc1)), ...
        min(tFirst(cc2))), CC(E(:,1)), CC(E(:,2)));
    
    tOverlapLast = cellfun(@(cc1,cc2) min(max(tLast(cc1)), ...
        max(tLast(cc2))), CC(E(:,1)), CC(E(:,2)));
    
    overlap = tOverlapLast - tOverlapFirst + 1;
        
    % Compute model parameters of each CC
    allCCPairParams1 = computeCCPairParams(allTrackParams,tFirst,tLast,...
        pFirst,CC(E(:,1)),tOverlapFirst,tOverlapLast,overlap);
    
    allCCPairParams2 = computeCCPairParams(allTrackParams,tFirst,tLast,...
        pFirst,CC(E(:,2)),tOverlapFirst,tOverlapLast,overlap);
    
    % Compute the edge weight
    allCCPairParams1 = num2cell(allCCPairParams1,1);
    [x1,y1,sx1,t1] = allCCPairParams1{:};
    
    allCCPairParams2 = num2cell(allCCPairParams2,1);
    [x2,y2,sx2,t2] = allCCPairParams2{:};

    % WD is the weight associated with pairwise distance between pair of CC
    allPWD = sqrt((x1-x2).^2 + (y1-y2).^2);
    WD = (allPWD + 1).^(-1/2);
    
    % WA is the weight associated with the deviation of each feature from
    % the pair axis
    ct = cos(t1);
    st = sin(t1);
    x0 = x1 - sx1 .* ct;
    y0 = y1 - sx1 .* st;
    pD1 = abs(((x2-x1) .* (y1 - y0) - (x1 - x0) .* (y2 - y1)) ./ allPWD);
    
    ct = cos(t2);
    st = sin(t2);
    x0 = x2 - sx2 .* ct;
    y0 = y2 - sx2 .* st;
    pD2 = abs(((x2-x1) .* (y1 - y0) - (x1 - x0) .* (y2 - y1)) ./ allPWD);
    
    WA = exp(-.5 * (pD1.^2 + pD2.^2));
    
    ppLast = cumsum(overlap);
    ppFirst = ppLast-overlap+1;
    
    W = arrayfun(@(a,b) mean(WA(a:a+b-1) .* WD(a:a+b-1)), ppFirst, overlap);
    
    % Compute the pairwise matching
    M = maxWeightedMatching(nCC,E,W);

    % There are 3 categories of CC pairs:
    % - unmatching pairs: keep them in E and create
    % - significant matched pairs: remove them from E and create new
    % - unsignificant matched pairs: remove them once and for all

    % Threshold
    % TODO
    
    % Update CC and other related variables
    % TODO
    
    % Save CC into a file
    % Save labeled tracks into a file
    % TODO
end

%% END
movieData.pairTracks.status = 1;
