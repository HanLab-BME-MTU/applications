function movieData = getMoviePairTracks(movieData,maxDistance,...
    minOverlap,bandWidth,minDistance,alpha,batchMode)

% maxDistance:     maximum distance between track pair candidates (t1,t2) such
%                  that max(dist(t1,t2)) < maxDistance.
%
% minOverlap:      is the minimal number of frame 2 tracks need to live together
%                  to be considered as potential pair of tracks. Default is 1.
%
% bandWidth:       distance in nanometers away from the cell edge where pair of
%                  feature needs to be evaluate whether it follows the Actin
%                  flow.
%
% minDistance:     minimum distance between 2 track pair candidates in the
%                  first 'bandWidth' nanometers away from cell edge.
%
% alpha:           quantile of PDF tail. Default is 0.05.

% Parse input parameters

if nargin < 3 || isempty(minOverlap)
    minOverlap = 2;
end

assert(minOverlap >= 2);

if nargin < 4 || isempty(bandWidth)
    bandWidth = 1000;
end

if nargin < 5 || isempty(alpha)
    alpha = 0.05;
end

assert(checkMovieDistanceTransform(movieData));
assert(checkMovieParticleDetection(movieData));
assert(checkMovieParticleTracking(movieData));

% BEGIN
movieData.pairTracks.status = 0;

% Create output directory
movieData.pairTracks.directory = fullfile(movieData.analysisDirectory, 'pairTracks');

if ~exist(movieData.pairTracks.directory, 'dir')
    mkdir(movieData.pairTracks.directory);
end

% Get all track parameters
[allTrackParams, tFirst, lifetime] = getAllTrackParams(movieData);

% Get initial pair track candidates
[E, overlap, pFirst1, pFirst2] = ...
    getInitialPairTracks(movieData, allTrackParams, tFirst, lifetime, ...
    maxDistance, minOverlap, bandWidth, minDistance);

% For each pair i = (t1, t2), we have:
% - trackModels1(i,:)    = segment model parameters of every points of
%                          track t1 over overlap(i) frames
% - trackModels2(i,:)    = segment model parameters of every points of
%                          track t2 over overlap(i) frames
% - pairTrackModels(i,:) = segment model parameters of every point of track
%                          t1 AND t2 over overlap(i) frames

% Compute initial model for each first track within pair candidates
allTrackParams1 = arrayfun(@(a,b) allTrackParams(a:a+b-1,:), ...
    pFirst1, overlap, 'UniformOutput', false);
allTrackParams1 = vertcat(allTrackParams1{:});
[trackModels1 res1] = getTrackModels(allTrackParams1, overlap);

% Compute initial model for each second track within pair candidates
allTrackParams2 = arrayfun(@(a,b) allTrackParams(a:a+b-1,:), ...
    pFirst2, overlap, 'UniformOutput', false);
allTrackParams2 = vertcat(allTrackParams2{:});
[trackModels2 res2] = getTrackModels(allTrackParams2, overlap);

% Compute initial model for each pair candidates (merged tracks)
pLast = cumsum(overlap);
pFirst = pLast-overlap+1;
allPairTrackParams = arrayfun(@(a,b) [allTrackParams1(a:a+b-1,:); ...
    allTrackParams2(a:a+b-1,:)], pFirst, overlap, 'UniformOutput', false);
allPairTrackParams = vertcat(allPairTrackParams{:});
[pairTrackModels resPair] = getTrackModels(allPairTrackParams, 2 * overlap);

% Compute BIC of the split model versus merged model
N = 4 * overlap;
varSplit = cellfun(@(res1,res2) var([res1; res2], 1), res1, res2);
varMerge = cellfun(@(res) var(res,1), resPair);

bicSplit = N .* log(varSplit) + 8 * log(N);
bicMerge = N .* log(varMerge) + 4 * log(N);

% Define score of each pair
W = bicSplit - bicMerge;

% Define which pair is valid
isValidPair = false(size(E,1),1);
isValidPair(W >= 0) = true;

% DEBUG
nFrames = movieData.nImages(1);
tOverlapFirst = max(tFirst(E(:,1)), tFirst(E(:,2)));

segments = cell(nFrames,1);
for iFrame = 1:nFrames
    isPairInFrame = iFrame >= tOverlapFirst & ...
        iFrame <= tOverlapFirst + overlap - 1;
    
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

save(fullfile(movieData.pairTracks.directory, 'allPairTrackCands.mat'), 'segments');

% DEBUG
segments = cell(nFrames,1);
for iFrame = 1:nFrames
    isPairInFrame = iFrame >= tOverlapFirst & ...
        iFrame <= tOverlapFirst + overlap - 1 & ...
        isValidPair;
    
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

save(fullfile(movieData.pairTracks.directory, 'allValidPairTrackCands.mat'), 'segments');


% % Start iteration
% while any(isValidPair)
%     M = maxWeightedMatching(E(isValidPair,:), W(isValidPair));
% end

% END
movieData.pairTracks.status = 1;
