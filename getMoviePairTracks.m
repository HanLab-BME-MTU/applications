function movieData = getMoviePairTracks(movieData,minOverlap,timeGap,maxEuclidianDist,sigmaPSF,nLevels,alpha,probBinSize,batchMode)

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
% sigmaPSF:             the standard deviation of a Gaussian model of PSF
%
% nLevels:              maximum number of pair level.
%
% alpha:                quantile of PDF tail. Default is 0.05.
%
% probBinSize:          size in [0...1] of a bin. probBinSize^-1 gives the
%                       the number of bins a PDF is cut into. Default is
%                       1e-4.

% 0) parse input parameters

if nargin < 2 || isempty(minOverlap)
    minOverlap = 1;
end

if nargin < 3 || isempty(timeGap)
    timeGap = 0;
end

if nargin < 6 || isempty(nLevels)
    nLevels = 5;
end

if nargin < 7 || isempty(alpha)
    alpha = 0.05;
end

if nargin < 8 || isempty(probBinSize)
    probBinSize = 1e-4;
end

assert(checkMovieParticleDetection(movieData));
assert(checkMovieParticleTracking(movieData));

movieData.pairTracks.status = 0;

movieData.pairTracks.directory = fullfile(movieData.analysisDirectory, 'pairTracks');

% Create output directory
if ~exist(movieData.pairTracks.directory, 'dir')
    mkdir(movieData.pairTracks.directory);
end

imagePath = fullfile(movieData.directory, movieData.channelDirectory{1});
imageFiles = dir([imagePath filesep '*.tif']);

nFrames = movieData.nImages(1);

% 1) Preprocess tracks

load(fullfile(movieData.particleDetection.directory, movieData.particleDetection.filename));
load(fullfile(movieData.particleTracking.directory, movieData.particleTracking.filename));

nTracks = size(tracksFinal,1); %#ok<NODEF>

% Check there is no split and merge
% if split and merge was enabled, it would complexify the interpolation in
% gaps (see step 3)
seqOfEvents = vertcat(tracksFinal(:).seqOfEvents);
assert(nnz(isnan(seqOfEvents(:,4))) == size(seqOfEvents,1));

SEL = getTrackSEL(tracksFinal);
iFirst = SEL(:,1)';
iLast = SEL(:,2)';
lifetime = SEL(:,3)';

% 2) Interpolate feature parameters in gaps

frameTracks = arrayfun(@(a,b) a:b, iFirst, iLast, 'UniformOutput', false);
frameTracks = horzcat(frameTracks{:});
tracksFeatIndx = [tracksFinal(:).tracksFeatIndxCG];
% params = parameter values of all tracks:
% x11 dx11 y11 dy11 A11 dA11 Sx11 dSx11 Sy11 dSy11 T11 dT11
% x12 dx12 ... (2nd point of the 1st track)
% ...
% x1n dx1n ... (nth point of the 1st track)
% x21 dx21 ... (1st point of the 2nd track)
% ...
params = nan(numel(tracksFeatIndx),12);

for iFrame = 1:nFrames
    ind = frameTracks == iFrame & tracksFeatIndx ~= 0;
    xCoord = [featuresInfo(iFrame).xCoord];
    yCoord = [featuresInfo(iFrame).yCoord];
    amp = [featuresInfo(iFrame).amp];
    sX = [featuresInfo(iFrame).stdAlong];
    sY = [featuresInfo(iFrame).stdAside];
    theta = [featuresInfo(iFrame).theta];
    
    params(ind,1:2) = xCoord(tracksFeatIndx(ind),:);
    params(ind,3:4) = yCoord(tracksFeatIndx(ind),:);
    params(ind,5:6) = amp(tracksFeatIndx(ind),:);
    params(ind,7:8) = sX(tracksFeatIndx(ind),:);
    params(ind,9:10) = sY(tracksFeatIndx(ind),:);
    params(ind,11:12) = theta(tracksFeatIndx(ind),:);
end

gacombIdx = diff(isnan(params(:,1)));
gapStarts = find(gacombIdx==1)+1;
gapEnds = find(gacombIdx==-1);
gapLengths = gapEnds-gapStarts+1;
nGaps = length(gapLengths);

for g = 1:nGaps
    borderIdx = [gapStarts(g)-1 gapEnds(g)+1];
    gacombIdx = gapStarts(g):gapEnds(g);
    params(gacombIdx,1:2) = interp1(borderIdx, params(borderIdx,1:2), gacombIdx);
    params(gacombIdx,3:4) = interp1(borderIdx, params(borderIdx,3:4), gacombIdx);
    params(gacombIdx,5:6) = interp1(borderIdx, params(borderIdx,5:6), gacombIdx);
    params(gacombIdx,7:8) = interp1(borderIdx, params(borderIdx,7:8), gacombIdx);
    params(gacombIdx,9:10) = interp1(borderIdx, params(borderIdx,9:10), gacombIdx);
    params(gacombIdx,11:12) = interp1(borderIdx, params(borderIdx,11:12), gacombIdx);
end

% 3) Compute the set of valid pairs of tracks that significantly overlap in
% time.

% First, end indexes for X and Y
last = cumsum(lifetime);
first = last-lifetime+1;

validTrackPairIdx = pcombs(1:nTracks);

overlapFirst = max(iFirst(validTrackPairIdx(:,1)), iFirst(validTrackPairIdx(:,2)));
overlapLast = min(iLast(validTrackPairIdx(:,1)), iLast(validTrackPairIdx(:,2)));
overlap = overlapLast - overlapFirst + 1 + timeGap;

hasOverlap = overlap >= minOverlap;

% trim arrays
validTrackPairIdx = validTrackPairIdx(hasOverlap,:);
overlapFirst = overlapFirst(hasOverlap);
overlapLast = overlapLast(hasOverlap);
overlap = overlap(hasOverlap);

% 4) Compute the euclidian distance between pair of tracks

% translate overlap first/last values to 1-D indexes
first1 = first(validTrackPairIdx(:,1)) + overlapFirst - iFirst(validTrackPairIdx(:,1));
first2 = first(validTrackPairIdx(:,2)) + overlapFirst - iFirst(validTrackPairIdx(:,2));

% sort overlap values
[overlap idx] = sort(overlap);
validTrackPairIdx = validTrackPairIdx(idx,:);

first1 = first1(idx);
first2 = first2(idx);

firstIdx = find([1 diff(overlap)]);
lastIdx = find([-diff(-overlap) 1]);

overlapValues = unique(overlap);

euclidianDist = zeros(size(overlap));

for k = 1:length(firstIdx)
    % indexes corresponding to overlap value
    range = firstIdx(k):lastIdx(k);
    
    M = repmat((0:overlapValues(k)-1), [length(range) 1]);
    idx1 = repmat(first1(range)', [1 overlapValues(k)]) + M;
    idx2 = repmat(first2(range)', [1 overlapValues(k)]) + M;
    
    x1 = params(idx1,1);
    y1 = params(idx1,3);
    x2 = params(idx2,1);
    y2 = params(idx2,3);
    
    % average distance
    euclidianDist(range) = mean(reshape(sqrt((x1-x2).^2 + (y1-y2).^2), ...
        [length(range) overlapValues(k)]), 2);
end

% 5) Trim the pair of tracks that are too far apart from each other
isCloseEnough = euclidianDist <= maxEuclidianDist;
validTrackPairIdx = validTrackPairIdx(isCloseEnough,:);

% Save pair tracks per frame
trackLabels = computeLabel(validTrackPairIdx, nTracks); %#ok<NASGU>
save(fullfile(movieData.pairTracks.directory, ...
    'classifiedTracksPerDistance.mat'),'tracksFinal', 'trackLabels');

% CC: a cell array of connected components. Each connected component
% contains a list of track indices
uniqueTrackIdx = unique(validTrackPairIdx(:));
CC = arrayfun(@(x) {x}, uniqueTrackIdx);
nCC = numel(CC);

% CCparams: a struct array of size = numel(CC). Each element is the
% model parameters representing a particular CC.
CCparams(1:nCC) = struct('iFirst',[],'iLast',[],'modelType',[],'params',[]);
% CCparams
%   .iFirst:     birth frame of the CC
%   .iLast:      death frame of the CC
%   .modelType:  0 = anisotropic gaussian model, 1 = segment model
%   .params:     model parameter organized as follows:
%      if modelType == 0:
%            x1 dx1 y1 dy1 A1 dA1 Sx1 dSx1 Sy1 dSy1 T1 dT1  (model
%            parameters at frame iFirst)
%            x2 dx2 y2 dy2 A2 dA2 Sx2 dSx2 Sy2 dSy2 T2 dT2  (model
%            parameters at frame iFirst+1)
%            ...
%            xN dxN yN ... (model parameters at frame iLast
%
%      if modelType == 1:
%            TBD
%   .bounds:     xmin,ymin,xmax,ymax
for iCC = 1:nCC
    iTrack = CC{iCC};
    CCparams(iCC).iFirst = iFirst(iTrack);
    CCparams(iCC).iLast = iLast(iTrack);
    CCparams(iCC).modelType = 0;
    CCparams(iCC).params = params(first(iTrack):last(iTrack),:);
end

% CCpairIdx: pair of CC indices for which we need a weight. We initialize
% CCpairIdx with validTrackPairIdx
tmp = zeros(max(uniqueTrackIdx),1);
tmp(uniqueTrackIdx) = 1:numel(uniqueTrackIdx);
CCpairIdx = reshape(tmp(validTrackPairIdx(:)), size(validTrackPairIdx));
nCCPairs = size(CCpairIdx,1);

for iLevel = 1:nLevels

    % Initialized model of each pair of CC    
    modelType1 = [CCparams(CCpairIdx(:,1)).modelType];
    modelType2 = [CCparams(CCpairIdx(:,2)).modelType];
    
    isCCpair00 = ~(modelType1 | modelType2);
    isCCpair01 = ~modelType1 & modelType2;
    isCCpair10 = modelType1 & ~modelType2;
    isCCpair11 = modelType1 & modelType2;

    CCpairParams(1:nCCPairs) = struct(...
        'iFirst',[],...
        'iLast',[],...
        'modelType',[],...
        'params',[],...
        'bounds',[]);
    
    CCpairParams(isCCpair00) = arrayfun(@(X,Y) initPairParams00(X,Y),...
        CCparams(CCpairIdx(isCCpair00,1)),...
        CCparams(CCpairIdx(isCCpair00,2)));
    
    CCpairParams(isCCpair01) = arrayfun(@(X,Y) initPairParams01(X,Y),...
        CCparams(CCpairIdx(isCCpair01,1)),...
        CCparams(CCpairIdx(isCCpair01,2)));
    
    CCpairParams(isCCpair10) = arrayfun(@(X,Y) initPairParams01(X,Y),...
        CCparams(CCpairIdx(isCCpair10,2)),...
        CCparams(CCpairIdx(isCCpair10,1)));

    CCpairParams(isCCpair11) = arrayfun(@(X,Y) initPairParams11(X,Y),...
        CCparams(CCpairIdx(isCCpair11,1)),...
        CCparams(CCpairIdx(isCCpair11,2)));
    
    % Fit model for each pair of CC
    iFirst = [CCpairParams(:).iFirst];
    iLast = [CCpairParams(:).iLast];

    for iFrame = 1:nFrames
        % Read image
        ima = imread(fullfile(imagePath, imageFiles(iFrame).name));
        
        % Get which pair of CC is living at frame iFrame
        indInFrame = find(iFirst <= iFrame & iLast >= iFrame);
        
        for iPair = indInFrame
        
            offset = iFrame - iFirst(iPair) + 1;
        
            initParams = CCpairParams(iPair).params(offset,:);
        
            bounds = CCpairParams(iPair).bounds;
            
            crop = ima(bounds(2):bounds(4), bounds(1):bounds(3));
            
            CCpairParams(iPair).params(offset,:) = ...
                fitSegment2D(crop, initParams, 'xyArltC'); % TODO
        end
    end
    
    % Test CC pairs candidate (over the overlapping period):
    % - extreme deviations from CCParams
    % - goodness-of-fit
    % - model parameters
    % TODO

    % TODO: remove any pair of CC from CCpairIdx for which the tests have
    % failed.
    % DO NOT FORGET TO UPDATE ANY ARRAY RELATED TO CCPAIRIDX !!!
    
    % Compute weights for each pair of CC
    
    W = zeros(nCCPairs,1);
    % TODO
    
    % Compute Maximum Matching
    D = sparse(CCpairIdx(:,1),CCpairIdx(:,2),W,nCC,nCC,numel(W));
    M = maxWeightMatching(D,eps);
    
    % Accept CCpairParams(M) and merge CCpairIdx(M,1) with CCpairIdx(M,2)
    % TODO
    
    % merge CCparams likewise and replace the merged CC params with
    % CCPairParams
    % TODO
    
    % Get a new set of CC pairs candidates
    % TODO
end

movieData.pairTracks.status = 1;

function pairParams = initPairParams00(params1,params2)

pairParams = struct(...
    'iFirst',[],...
    'iLast',[],...
    'modelType',[],...
    'params',[],...
    'bounds',[]);

pairParams.iFirst = max(params1.iFirst, params2.iFirst);
pairParams.iLast = min(params1.iLast, params1.iLast);
pairParams.modelType = 1;
% TODO !!!

function pairParams = initPairParams01(params1,params2)

pairParams = struct(...
    'iFirst',[],...
    'iLast',[],...
    'modelType',[],...
    'params',[],...
    'bounds',[]);

pairParams.iFirst = max(params1.iFirst, params2.iFirst);
pairParams.iLast = min(params1.iLast, params1.iLast);
pairParams.modelType = 1;
% TODO !!!

function pairParams = initPairParams11(params1,params2)

pairParams = struct(...
    'iFirst',[],...
    'iLast',[],...
    'modelType',[],...
    'params',[],...
    'bounds',[]);

pairParams.iFirst = max(params1.iFirst, params2.iFirst);
pairParams.iLast = min(params1.iLast, params1.iLast);
pairParams.modelType = 1;
% TODO !!!

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

