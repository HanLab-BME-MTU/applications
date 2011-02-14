function movieData = getMoviePairTracks(movieData,minOverlap,timeGap,maxEuclidianDist,sigmaPSF,accT,alpha,probBinSize,batchMode)

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
% accT:                 angle accuracy as a portion of pi. Default is 1/16.
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

if nargin < 6 || isempty(accT)
    accT = 1/16;
end

if nargin < 7 || isempty(alpha)
    alpha = 0.05;
end

if nargin < 8 || isempty(probBinSize)
    probBinSize = 1e-4;
end

assert(checkMovieParticleTracking(movieData));

movieData.pairTracks.status = 0;

movieData.pairTracks.directory = fullfile(movieData.analysisDirectory, 'pairTracks');

% Create output directory
if ~exist(movieData.pairTracks.directory, 'dir')
    mkdir(movieData.pairTracks.directory);
end

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
X = nan(numel(tracksFeatIndx),2);
Y = nan(numel(tracksFeatIndx),2);
A = nan(numel(tracksFeatIndx),2);
Sx = nan(numel(tracksFeatIndx),2);
Sy = nan(numel(tracksFeatIndx),2);
T = nan(numel(tracksFeatIndx),2);

for iFrame = 1:nFrames
    ind = frameTracks == iFrame & tracksFeatIndx ~= 0;
    xCoord = [featuresInfo(iFrame).xCoord];
    yCoord = [featuresInfo(iFrame).yCoord];
    amp = [featuresInfo(iFrame).amp];
    sX = [featuresInfo(iFrame).stdAlong];
    sY = [featuresInfo(iFrame).stdAside];
    theta = [featuresInfo(iFrame).theta];
    
    X(ind,:) = xCoord(tracksFeatIndx(ind),:);
    Y(ind,:) = yCoord(tracksFeatIndx(ind),:);
    A(ind,:) = amp(tracksFeatIndx(ind),:);
    Sx(ind,:) = sX(tracksFeatIndx(ind),:);
    Sy(ind,:) = sY(tracksFeatIndx(ind),:);
    T(ind,:) = theta(tracksFeatIndx(ind),:);
end

gacombIdx = diff(isnan(X(:,1)));
gapStarts = find(gacombIdx==1)+1;
gapEnds = find(gacombIdx==-1);
gapLengths = gapEnds-gapStarts+1;
nGaps = length(gapLengths);

for g = 1:nGaps
    borderIdx = [gapStarts(g)-1 gapEnds(g)+1];
    gacombIdx = gapStarts(g):gapEnds(g);
    X(gacombIdx,:) = interp1(borderIdx, X(borderIdx,:), gacombIdx);
    Y(gacombIdx,:) = interp1(borderIdx, Y(borderIdx,:), gacombIdx);
    A(gacombIdx,:) = interp1(borderIdx, A(borderIdx,:), gacombIdx);
    Sx(gacombIdx,:) = interp1(borderIdx, Sx(borderIdx,:), gacombIdx);
    Sy(gacombIdx,:) = interp1(borderIdx, Sy(borderIdx,:), gacombIdx);
    T(gacombIdx,:) = interp1(borderIdx, T(borderIdx,:), gacombIdx);
end

% 3) Find the set of track pair candidates that significantly overlap in
% time.

% First, end indexes for X and Y
last = cumsum(lifetime);
first = last-lifetime+1;

pairIdx = pcombs(1:nTracks);

overlapFirst = max(iFirst(pairIdx(:,1)), iFirst(pairIdx(:,2)));
overlapLast = min(iLast(pairIdx(:,1)), iLast(pairIdx(:,2)));
overlap = overlapLast - overlapFirst + 1 + timeGap;

hasOverlap = overlap >= minOverlap;

fprintf('Overlapping track pairs = %f %%\n',...
    nnz(hasOverlap) * 100 / numel(hasOverlap));

% trim arrays
pairIdx = pairIdx(hasOverlap,:);
overlapFirst = overlapFirst(hasOverlap);
overlapLast = overlapLast(hasOverlap); %#ok<NASGU>
overlap = overlap(hasOverlap);

% 4) Compute the euclidian distance between pair of tracks

% translate overlap first/last values to 1-D indexes
first1 = first(pairIdx(:,1)) + overlapFirst - iFirst(pairIdx(:,1));
first2 = first(pairIdx(:,2)) + overlapFirst - iFirst(pairIdx(:,2));

% sort overlap values
[overlap idx] = sort(overlap);
pairIdx = pairIdx(idx,:);

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
    
    x1 = X(idx1);
    y1 = Y(idx1);
    x2 = X(idx2);
    y2 = Y(idx2);
    
    % average distance
    euclidianDist(range) = mean(reshape(sqrt((x1-x2).^2 + (y1-y2).^2), ...
        [length(range) overlapValues(k)]), 2);
end

% 5) Trim the pair of tracks that are too far apart from each other
isCloseEnough = euclidianDist <= maxEuclidianDist;
fprintf('Close track pairs = %f %%\n',...
    nnz(isCloseEnough) * 100 / numel(hasOverlap));
pairIdx = pairIdx(isCloseEnough,:);

% Save pair tracks per frame
trackLabels = computeLabel(pairIdx, nTracks); %#ok<NASGU>
save(fullfile(movieData.pairTracks.directory, ...
    'classifiedTracksPerDistance.mat'),'tracksFinal', 'trackLabels');

% CC: a cell array of connected components. Each connected component
% contains a list of track indices
CC = ;
% CCparams: a struct array of size = numel(CC). Each element is the
% model parameters representing a particular CC.
CCparams = ;
% CCpairIdx: pair of CC indices for which we need a weight
CCpairIdx = ;
% allowedTraiPairIdx: the list of allowed pair between tracks. Over the
% course of the algorithm, this list is decreasing
allowedTrackPairIdx = CCpairIdx;

for iLevel = 1:nLevels
    % Compute the weights between pair of CC
    [W, CCPairParams, allowedTrackPairIdx] = ...
        computeWeights(CC, CCparams, CCpairIdx, allowedTrackPairIdx);
    
    % Maximum Matching
    M = maxWeightMatching(CCpairIdx,W);
    
    % merge CCpairIdx(M,1) with CCpairIdx(M,2)
    
    % merge CCparams likewise and replace the mergegCC params with
    % CCPairParams
    
    % merge CC
    
    CCpairIdx = computePairCC(CCpairIdx, allowedTrackPairIdx);
end




% 9) Solve the matching problem

weights(weights == 0) = eps;

D = sparse(pairIdx(:,1),pairIdx(:,2), weights, nTracks, nTracks, size(pairIdx,1));

M = maxWeightMatching(D, eps);

% Save pair tracks per frame
trackLabels = computeLabel(M, nTracks); %#ok<NASGU>
save(fullfile(movieData.pairTracks.directory, ...
    'pairTrack.mat'),'tracksFinal', 'trackLabels');

movieData.pairTracks.status = 1;

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
% end
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

