function E = getInitialPairTracks1(nFrames,imSize,pixelSize, distToEdgePath, allFeatures, tFirst, ...
    lifetime, maxDistance, minOverlap, bandWidth, minDistance)

%This function requires the distance transform 

maxDistance = maxDistance / pixelSize;
bandWidth   = bandWidth / pixelSize;
minDistance = minDistance / pixelSize;

nTracks = numel(tFirst);
tLast   = tFirst + lifetime - 1;

% Discard any pair (t1,t2) such that
% max(dist(t1,t2)) >= maxDistance
% OR
% overlap(t1,t2) < minOverlap

% pFirst and pLast are indexing every variable named 'allTrack*'
pLast  = cumsum(lifetime);
pFirst = pLast-lifetime+1;

D = sparse(nTracks,nTracks);
N = sparse(nTracks,nTracks);

for iFrame = 1:nFrames
    indTrackInFrame = find(iFrame >= tFirst & iFrame < tFirst + lifetime);
    
    indP = pFirst(indTrackInFrame) + iFrame - tFirst(indTrackInFrame);
    X = allFeatures(indP, 1);
    Y = allFeatures(indP, 2);

    [indTrack1 dist] = KDTreeBallQuery([X,Y],[X,Y], repmat(maxDistance, numel(X),1));
    
    % remove self references (d == 0 which is always the first element)
    indTrack1 = cellfun(@(c) indTrackInFrame(c(2:end)), ...
        indTrack1, 'UniformOutput', false);
    
    dist = cellfun(@(c) c(2:end), dist, 'UniformOutput', false);
    dist = vertcat(dist{:});
    
    indTrack2 = arrayfun(@(a,b) repmat(a,b,1), indTrackInFrame, ...
        cellfun(@numel,indTrack1), 'UniformOutput', false);

    pairs = [vertcat(indTrack1{:}), vertcat(indTrack2{:})];
    
    % keep only those where pair(i,1) < pair(i,2)
    isValid  = pairs(:,1) < pairs(:,2);
    indPairs = sub2ind(size(D), pairs(isValid,1), pairs(isValid,2));
    dist     = dist(isValid);

    % FIXME !!!
    D(indPairs) = max(D(indPairs), dist);
    N(indPairs) = N(indPairs) + 1;
end

isValid = D ~= 0;
[I J]   = find(isValid);
E       = [I J];
[E idx] = sortrows(E);

D = full(D(isValid));
D = D(idx);
N = full(N(isValid));
N = N(idx);

% Compute the overlap between pair of tracks
tOverlapFirst = max(tFirst(E(:,1)), tFirst(E(:,2)));
tOverlapLast  = min(tLast(E(:,1)), tLast(E(:,2)));
overlap       = tOverlapLast - tOverlapFirst + 1;

% overlap == N means that over the overlap, (t1,t2) was always < maxDistance
isValid = overlap == N & overlap >= minOverlap;

% trim arrays
E = E(isValid,:);
tOverlapFirst = tOverlapFirst(isValid);
D = D(isValid);

% Point the location of each track parameters for every pair
pFirst1 = pFirst(E(:,1)) + tOverlapFirst - tFirst(E(:,1));
pFirst2 = pFirst(E(:,2)) + tOverlapFirst - tFirst(E(:,2));

% DEBUG
% filename = fullfile(movieData.pairTracks.directory, 'allPairTrackCands.mat');
% savePairTrack(movieData, allFeatures, tFirst, E, pFirst1, pFirst2, filename);

% Discard any pair (t1,t2) such that:
% distToEdge(t1) < bandWidth & distToEdge(t2) < bandWidth
% AND
% max(dist(t1,t2)) < minDistance


distToEdgeFiles = dir([distToEdgePath filesep '*.mat']);

isTrackInBand = false(nTracks,1);

for iFrame = 1:nFrames
    load(fullfile(distToEdgePath, distToEdgeFiles(iFrame).name));
    
    isTrackInFrame = iFrame >= tFirst & iFrame <= tLast;
    
    indP = pFirst(isTrackInFrame) + iFrame - tFirst(isTrackInFrame);
    X    = max(1,min(round(allFeatures(indP, 1)),imSize(1)));
    Y    = max(1,min(round(allFeatures(indP, 2)),imSize(2)));
    ind  = sub2ind(imSize,Y, X);
        
    isTrackInBand(isTrackInFrame) = isTrackInBand(isTrackInFrame) | ...
        distToEdge(ind) < bandWidth;
end

isTrackPairInBand = isTrackInBand(E(:,1)) & isTrackInBand(E(:,2));
isValid           = ~(isTrackPairInBand & D < minDistance);

% trim arrays
E       = E(isValid,:);
pFirst1 = pFirst1(isValid);
pFirst2 = pFirst2(isValid);

% DEBUG
% filename = fullfile(movieData.pairTracks.directory, 'pairTrackCandsInBands.mat');
% savePairTrack(movieData, allFeatures, tFirst, E, pFirst1, pFirst2, filename);
