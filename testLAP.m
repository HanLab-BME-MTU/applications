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

% Load images
nFrames = movieData.nImages;
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

% DEBUG: save tracks classification
fprintf('Feature-based track classification: static = %f %%, linear = %f %%\n',...
    nnz(~isFeatureAligned) * 100 / numel(isFeatureAligned), ...
    nnz(isFeatureAligned) * 100 / numel(isFeatureAligned));

trackLabels = isFeatureAligned + 1; %#ok<NASGU> % (1 == static, 2 == aligned)
save(fullfile(movieData.particleTracking.directory, ...
    'classifiedTracksPerFeatureAlignement.mat'),'tracksFinal', 'trackLabels');

% 4) Classify each track whether its dynamics linear or not

% thetaTracks is a cell array of size nTracks, where each cell contains an
% array of size = lifetime of the track -1 and contains the direction from
% one track point to the next
thetaTracks = arrayfun(@(track) zeros(1, size(track.tracksFeatIndxCG,2) - 1), ...
    tracksFinal, 'UniformOutput', false);

for iFrame = 1:nFrames
    % TODO
end

isDynamicsAligned = false(nTracks,1);

for iTrack = 1:nTracks
    theta = thetaTracks{iTrack};
    
    nonNaNTheta = theta(~isnan(theta));
    
    isInWedge = bsxfun(@(t,w) t >= w & t < w + pi * accT, nonNaNTheta', wedges);
    
    k = max(sum(isInWedge,1));

    isDynamicsAligned(iTrack) = k > cutoffs(numel(theta));
end

% DEBUG: save tracks classification
fprintf('Dynamics-based track classification: static = %f %%, linear = %f %%\n',...
    nnz(~isDynamicsAligned) * 100 / numel(isDynamicsAligned), ...
    nnz(isDynamicsAligned) * 100 / numel(isDynamicsAligned));

trackLabels = isDynamicsAligned + 1; %#ok<NASGU> % (1 == static, 2 == aligned)
save(fullfile(movieData.particleTracking.directory, ...
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

% DEBUG: save pair tracks per frame
trackLabels = computeLabel(pairIdx, nTracks); %#ok<NASGU>
save(fullfile(movieData.particleTracking.directory, ...
    'classifiedTracksPerDistance.mat'),'tracksFinal', 'trackLabels');

% 8) For every possible combination of pair type, compute the probability
% of association:
%
% STATIC / STATIC
% LINEAR / STATIC
% LINEAR / LINEAR

% 9) Trim probability < 1%

% 10) Display


% 3) Trim the set of pair candidates by assessing how far they are from
% each other (radon distance)
%
% Q: What the threshold values (t and alpha)
%
% 4) Compute the max-weight matching problem on radon distance-based
% similarity function
%
% Q: what is the track-track similarity function?
% Q: is the double -> int quantification works?
%
% 5) Post-processing: remove pairs that are unsignificant
%
% Q: What is unsignificant?
%
% 6) Compute the similarity function for track-segment pair candidate and
% segment-segment pair candidate?
%
% Q: What is the track-segment similarity function?
%
% 7) Redo steps 4-6 until convergence
%
% Q: What is the stop criteria?


% thE: [0, +inf)
% thA: [0, pi]
% thP: [0, 1]

% imagePath = fullfile(movieData.imageDirectory, movieData.channelDirectory{1});
% imageFiles = dir([imagePath filesep '*.tif']);
% ima = imread(fullfile(imagePath, imageFiles(1).name));
% 
% load(fullfile(movieData.particleDetection.directory, ...
%     movieData.particleDetection.filename));
% 
% X = [featuresInfo(1).xCoord, featuresInfo(1).yCoord];
% ind = sub2ind(size(ima),X(:,2),X(:,1));
% N = size(ind,1);
% 
% [~, T] = steerableFiltering(double(ima),2,2);
% 
% Y = [cos(T(ind)), sin(T(ind))];
% 
% pair = pcombs(1:N,false);
% 
% u0 = X(pair(:,1),:);
% u1 = u0 + Y(pair(:,1),:);
% v0 = X(pair(:,2),:);
% v1 = v0 + Y(pair(:,2),:);
% 
% isValid = true(size(pair,1),1);
% 
% % euclidian distance [0...+inf]
% dE = sqrt(sum((u0 - v0).^2,2));
% isValid = isValid & dE <= thE;
% 
% % angle between u and v
% % dot = abs(sum((u1 - u0) .* (v1 - v0),2));
% % dot(dot > 1) = 1;
% % dot(dot < -1) = -1;
% % dA = acos(dot);
% % isValid = isValid & dA <= thA;
% %
% % mean distance of u1 and v1 projected on the line (u0,v0)
% dP1 = sqrt(1 - (sum((u1-u0) .* (v0-u0),2) ./ dE).^2);
% dP2 = sqrt(1 - (sum((v1-v0) .* (v0-u0),2) ./ dE).^2);
% dP = dP1 .* dP2;
% isValid = isValid & dP <= thP;
% %
% % cost = exp(- (dE .* (1/pi) .* dA .* dP));
% 
% cost = 1./ dE;
% 
% % Build cost matrix
% i = pair(isValid,1);
% j = pair(isValid,2);
% c = cost(isValid);
% 
% % Populate the lower triangular part only
% D = sparse(j, i, c, N, N, numel(c));
% 
% % Compute Maximum Weight Matching
% M = maxWeightMatching(D);
% 
% % Display result
% imshow(ima,[]);
% hold on;
% 
% B = zeros(N);
% ind = sub2ind([N N],M(:,1),M(:,2));
% B(ind) = 1;
% 
% line(X(:,1),X(:,2),'LineStyle','none', 'Marker', '.', 'Color', 'g');
% gplot(B,X,'r');
% quiver(X(:,1),X(:,2),Y(:,1),Y(:,2),0,'b');


