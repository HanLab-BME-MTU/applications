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
[allFeatures, tFirst, lifetime] = getAllFeatures(movieData);

% Get initial pair track candidates
[E, overlap, pFirst1, pFirst2] = ...
    getInitialPairTracks(movieData, allFeatures, tFirst, lifetime, ...
    maxDistance, minOverlap, bandWidth, minDistance);

% Define the set of connected components
CC = arrayfun(@(x) {x}, (1:nTracks)');
nCC = numel(CC);

% For each pair i = (t1, t2), we have:
% - ccModels1(i,:)    = segment model parameters of every points of track
%                       t1 over overlap(i) frames
% - ccModels2(i,:)    = segment model parameters of every points of track
%                       t2 over overlap(i) frames
% - pairCCModels(i,:) = segment model parameters of every point of track
%                       t1 AND t2 over overlap(i) frames

% Compute initial model for each first track within pair candidates
allFeaturesCC1 = arrayfun(@(a,b) allFeatures(a:a+b-1,[1 2 4 6]), ...
    pFirst1, overlap, 'UniformOutput', false);
allFeaturesCC1 = vertcat(allFeaturesCC1{:});
[modelsCC1 res1] = getSegmentModels(allFeaturesCC1, overlap);

% Compute initial model for each second track within pair candidates
allFeaturesCC2 = arrayfun(@(a,b) allFeatures(a:a+b-1,[1 2 4 6]), ...
    pFirst2, overlap, 'UniformOutput', false);
allFeaturesCC2 = vertcat(allFeaturesCC2{:});
[modelsCC2 res2] = getSegmentModels(allFeaturesCC2, overlap);

% Compute initial model for each pair candidates (merged tracks)
pLast = cumsum(overlap);
pFirst = pLast-overlap+1;
allPairTrackParams = arrayfun(@(a,b) [allFeaturesCC1(a:a+b-1,:); ...
    allFeaturesCC2(a:a+b-1,:)], pFirst, overlap, 'UniformOutput', false);
allPairTrackParams = vertcat(allPairTrackParams{:});
[pairCCModels resPair] = getSegmentModels(allPairTrackParams, 2 * overlap);

% Compute BIC of the split model versus merged model
N = 4 * overlap;
varSplit = cellfun(@(res1,res2) var([res1; res2], 1), res1, res2);
varMerge = cellfun(@(res) var(res,1), resPair);

bicSplit = N .* log(varSplit) + 8 * log(N);
bicMerge = N .* log(varMerge) + 4 * log(N);

% Define score of each pair
W = bicSplit - bicMerge;

% Define which pair is valid
isValidPair = W >= 0;

% DEBUG
% filename = fullfile(movieData.pairTracks.directory, 'allPairTrackCands.mat');
% savePairTrack(movieData, allFeatures, tFirst, E, pFirst1, pFirst2, filename);

E = E(isValidPair,:);
W = W(isValidPair);
pFirst1 = pFirst1(isValidPair);
pFirst2 = pFirst2(isValidPair);

% DEBUG
% filename = fullfile(movieData.pairTracks.directory, 'allValidPairTrackCands.mat');
% savePairTrack(movieData, allFeatures, tFirst, E, pFirst1, pFirst2, filename);

% Start iteration
while size(E,1)
    M = maxWeightedMatching(E, W);
        
    % Update E, CC and nCC
    EM = E(M,:);
    ccInd = 1:nCC;
    ccInd(EM(:,2)) = EM(:,1);
    ccIndUnique = unique(ccInd);
    values = 1:max(ccIndUnique);
    values(ccIndUnique) = 1:numel(ccIndUnique);
    ccInd = values(ccInd);

    E = ccInd(E);
    E = unique(E,'rows');
    E = E(E(:,1) ~= E(:,2),:);
    
    for iE = 1:size(EM,1)
        CC{EM(iE,1)} = horzcat(CC{EM(iE,1)}, CC{EM(iE,2)});
        CC{EM(iE,2)} = [];
    end
    
    isEmpty = cellfun(@isempty,CC);
    CC = CC(~isEmpty);
    nCC = numel(CC);

    % Compute the lifetime of each CC
    tFirstCC = cellfun(@(trackIdx) min(tFirst(trackIdx)), CC); % first frame of CC
    tLastCC = cellfun(@(trackIdx) max(tLast(trackIdx)), CC);   % last frame of CC
    lifetimeCC = tLastCC - tFirstCC + 1;                       % lifetime of CC
    
    % Compute the overlap between CC
    tOverlapFirst = max(tFirst(E(:,1)), tFirst(E(:,2)));
    tOverlapLast = min(tLast(E(:,1)), tLast(E(:,2)));
    overlap = tOverlapLast - tOverlapFirst + 1;

    % For each track in the first/second CC, compute the first feature
    % index that appears in the overlap
    pFirstCC1 = ;
    pFirstCC2 = ;
    
    % For each track in the first/second CC, compute the last feature
    % index that appears in the overlap
    pLastCC1 = ;
    pLastCC2 = ;
    
    % Gather every feature of each track in the first/second CC between
    % tOverlapFirst and tOverlapLast
    allFeaturesCC1 = cellfun(@(aa,bb) arrayfun(@(a,b) ...
        allFeatures(a:a+b-1, [1 2 4 6]), aa, bb, 'UniformOutput', false), ...
        ppFirst1, ppLast1, 'UniformOutput', false);
    
    allFeaturesCC1 = cellfun(@(c) vertcat(c{:}), allFeaturesCC1);
    
    allFeaturesCC2 = cellfun(@(aa,bb) arrayfun(@(a,b) ...
        allFeatures(a:a+b-1, [1 2 4 6]), aa, bb, 'UniformOutput', false), ...
        ppFirst2, ppLast2, 'UniformOutput', false);
    
    allFeaturesCC2 = cellfun(@(c) vertcat(c{:}), allFeaturesCC2);

    % Compute the number of features per each first/second CC
    nFeatures1 = cellfun(@(a) size(a,1), allFeaturesCC1);    
    allFeaturesCC1 = vertcat(allFeaturesCC1{:});
    
    nFeatures2 = cellfun(@(a) size(a,1), allFeaturesCC2);
    allFeaturesCC2 = vertcat(allFeaturesCC2{:});
   
    % Compute model of the each first/second CC
    [modelsCC1 res1] = getSegmentModels(allFeaturesCC1, nFeatures1);
    [modelsCC2 res2] = getSegmentModels(allFeaturesCC2, nFeatures2);
    
    
    
    % Recompute W
    
    % Threshold E, W, pFirst1, pFirst2
end

% END
movieData.pairTracks.status = 1;
