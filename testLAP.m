function testLAP(movieData,iChannel,minOverlap,timeGap,maxEuclidianDist,sigmaPSF,accT,alphaT)

% iChannel:             index of the filament channel
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
% sigmaPSF:             the standard deviation of a Gaussian model of PSF
%
% accT:                 angle accuracy as a portion of pi. Default is 1/16.
%
% alphaT:               quantile of the binomial tail. Default is 0.05.

% 0) parse input parameters

if nargin < 3 || isempty(minOverlap)
    minOverlap = 1;
end

if nargin < 4 || isempty(timeGap)
    timeGap = 0;
end

if nargin < 7 || isempty(accT)
    accT = 1/16;
end

if nargin < 8 || isempty(alphaT)
    alphaT = 0.05;
end

assert(checkMovieParticleTracking(movieData));

movieData.pairTracks.status = 0;

movieData.pairTracks.directory = fullfile(movieData.analysisDirectory, 'pairTracks');

% Create output directory
if ~exist(movieData.pairTracks.directory, 'dir')
    mkdir(movieData.pairTracks.directory);
end

nFrames = movieData.nImages(iChannel);

% Load images
imagePath = fullfile(movieData.imageDirectory, movieData.channelDirectory{iChannel});
imageFiles = dir([imagePath filesep '*.tif']);

% 1) Preprocess tracks

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

% 2) Interpolate position in gaps

X = arrayfun(@(t) t.tracksCoordAmpCG(1:8:end)', tracksFinal, 'UniformOutput', false);
X = vertcat(X{:});
Y = arrayfun(@(t) t.tracksCoordAmpCG(2:8:end)', tracksFinal, 'UniformOutput', false);
Y = vertcat(Y{:});

gacombIdx = diff(isnan(X));
gapStarts = find(gacombIdx==1)+1;
gapEnds = find(gacombIdx==-1);
gapLengths = gapEnds-gapStarts+1;
nGaps = length(gapLengths);

for g = 1:nGaps
    borderIdx = [gapStarts(g)-1 gapEnds(g)+1];
    gacombIdx = gapStarts(g):gapEnds(g);
    X(gacombIdx) = interp1(borderIdx, X(borderIdx), gacombIdx);
    Y(gacombIdx) = interp1(borderIdx, Y(borderIdx), gacombIdx);
end

% 3) Classify each track whether the direction given by the image feature
% (filtering) is constant along the track

% First, end indexes for X and Y
last = cumsum(lifetime);
first = last-lifetime+1;

% thetaTracks is a cell array of size nTracks, where each cell contains an
% array of size = lifetime of the track and contains the local orientaion
% along the track from each frame
thetaTracks = arrayfun(@(track) zeros(1, size(track.tracksFeatIndxCG,2)), ...
    tracksFinal, 'UniformOutput', false);
    
for iFrame = 1:nFrames
    % Get all the tracks that are in iFrame
    isInFrame = SEL(:,1) <= iFrame & SEL(:,2) >= iFrame;
    indTrack = find(isInFrame);
    
    indFeature = iFrame-SEL(isInFrame,1) + 1;
    
    % Get all positions
    x = X(first(indTrack) + indFeature' - 1);
    y = Y(first(indTrack) + indFeature' - 1);
        
    % Filter image
    ima = double(imread(fullfile(imagePath, imageFiles(iFrame).name)));
    [~,T] = steerableFiltering(ima,2,sigmaPSF);
    ct = cos(T);
    st = sin(T);
    
    % Interpolate vector component at (x,y) positions
    dU = interp2(ct, x, y,'cubic');
    dV = interp2(st, x, y,'cubic');
    
    % Normalize vector
    norm = sqrt(dU.^2 + dV.^2);     % norm cannot be 0
    dU = bsxfun(@rdivide,dU,norm);
    dV = bsxfun(@rdivide,dV,norm);
    
    theta = atan2(dU,dV);
    
    % Store theta
    for iiTrack = 1:numel(indTrack)
        iTrack = indTrack(iiTrack);
        thetaTracks{iTrack}(indFeature(iiTrack)) = theta(iiTrack);
    end
end

% compute the cutoff values so that the binomial tail = alphaTh
cutoffs = icdf('bino', 1-alphaT, 1:nFrames, accT);

isFeatureAligned = false(nTracks,1);

wedges = linspace(-pi/2, pi/2 - pi * accT, accT^-1 + 1);

for iTrack = 1:nTracks
    theta = thetaTracks{iTrack};
    
    isInWedge = bsxfun(@(t,w) t >= w & t < w + pi * accT, theta', wedges);
    
    k = max(sum(isInWedge,1));

    isFeatureAligned(iTrack) = k > cutoffs(numel(theta));
end

% Save tracks classification
percent = nnz(~isFeatureAligned) * 100 / numel(isFeatureAligned);
fprintf('Feature-based track classification: static = %f %%, linear = %f %%\n',...
    percent, 100 - percent);

trackLabels = isFeatureAligned + 1; %#ok<NASGU> % (1 == static, 2 == aligned)
save(fullfile(movieData.pairTracks.directory, ...
    'classifiedTracksPerFeatureAlignement.mat'),'tracksFinal', 'trackLabels');

% 4) Classify each track whether its dynamics linear or not

% thetaTracks is a cell array of size nTracks, where each cell contains an
% array of size = lifetime of the track -1 and contains the direction from
% one track point to the next
thetaTracks = cell(nTracks,1);

for iTrack = 1:nTracks
    dU = X(first(iTrack)+1:last(iTrack)) - X(first(iTrack):last(iTrack)-1);
    dV = Y(first(iTrack)+1:last(iTrack)) - Y(first(iTrack):last(iTrack)-1);
    
    % Normalize vector
    norm = sqrt(dU.^2 + dV.^2);     % norm can be 0 !!!
    
    theta = zeros(size(norm));
    isNull = abs(norm) < eps;
    theta(isNull) = NaN;
    
    dU = bsxfun(@rdivide,dU(~isNull),norm(~isNull));
    dV = bsxfun(@rdivide,dV(~isNull),norm(~isNull));
    theta(~isNull) = atan2(dU,dV);
    
    thetaTracks{iTrack} = theta;
end

isDynamicsAligned = false(nTracks,1);

for iTrack = 1:nTracks
    theta = thetaTracks{iTrack};
    
    nonNaNTheta = theta(~isnan(theta));

    if ~isempty(nonNaNTheta)
        isInWedge = bsxfun(@(t,w) t >= w & t < w + pi * accT, nonNaNTheta, wedges);
    
        k = max(sum(isInWedge,1));

        isDynamicsAligned(iTrack) = k > cutoffs(numel(theta));
    end
end

% Save tracks classification
percent = nnz(~isDynamicsAligned) * 100 / numel(isDynamicsAligned);
fprintf('Dynamics-based track classification: static = %f %%, linear = %f %%\n',...
    percent, 100 - percent);

trackLabels = isDynamicsAligned + 1; %#ok<NASGU> % (1 == static, 2 == aligned)
save(fullfile(movieData.pairTracks.directory, ...
    'classifiedTracksPerDynamicsAlignement.mat'),'tracksFinal', 'trackLabels');

% 5) Find the set of track pair candidates that significantly overlap in
% time.

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

% 6) Compute the euclidian distance between pair of tracks

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

% 7) Trim the pair of tracks that are too far apart from each other
pairIdx = pairIdx(euclidianDist <= maxEuclidianDist,:);

% Save pair tracks per frame
trackLabels = computeLabel(pairIdx, nTracks); %#ok<NASGU>
save(fullfile(movieData.pairTracks.directory, ...
    'classifiedTracksPerDistance.mat'),'tracksFinal', 'trackLabels');

% 8) For every possible combination of pair type, compute the probability
% of association:
%
% STATIC / STATIC => NO LINK !!!
% LINEAR / STATIC
% LINEAR / LINEAR

% 9) Trim probability < 1%

% 10) Display

movieData.pairTracks.status = 1;

