function testLAP(movieData,minOverlap,maxEuclidianDist)
% 1) Preprocess tracks

load(fullfile(movieData.particleTracking.directory, ...
    movieData.particleTracking.filename));

nTracks = size(tracksFinal,1); %#ok<NODEF>

% Check there is no split and merge
% if split and merge was enabled, it would complexify the interpolation in
% gaps (see step 3)
seqOfEvents = vertcat(tracksFinal(:).seqOfEvents);
assert(nnz(isnan(seqOfEvents(:,4))) == size(seqOfEvents,1));

% 2) Find the set of track pair candidates that significantly overlap in
% time.

SEL = getTrackSEL(tracksFinal);
iFirst = SEL(:,1)';
iLast = SEL(:,2)';
lifetime = SEL(:,3)';
pairIdx = pcombs(1:nTracks);

overlapFirst = max(iFirst(pairIdx(:,1)), iFirst(pairIdx(:,2)));
overlapLast = min(iLast(pairIdx(:,1)), iLast(pairIdx(:,2)));
overlap = overlapLast - overlapFirst + 1;

hasOverlap = overlap >= minOverlap;

fprintf('Overlapping track pairs = %f %%\n',...
    nnz(hasOverlap) * 100 / numel(hasOverlap));

% trim arrays
pairIdx = pairIdx(hasOverlap,:);
overlapFirst = overlapFirst(hasOverlap);
overlapLast = overlapLast(hasOverlap); %#ok<NASGU>
overlap = overlap(hasOverlap);

% 3) Interpolate position in gaps

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

% 4) Compute the euclidian distance between pair of tracks

% First, end indexes for X and Y
last = cumsum(lifetime);
first = last-lifetime+1;

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
pairIdx = pairIdx(euclidianDist <= maxEuclidianDist,:);

fprintf('Neighboring track pairs = %f %%\n',...
    size(pairIdx,1) * 100 / numel(hasOverlap));



% DEBUG: save pair tracks per frame


% Calculate the connected component and store the CC and the tracks
% (instead of storing the edges only).

G = sparse(pairIdx(:,1),pairIdx(:,2),true(size(pairIdx,1),1), ...
    nTracks, nTracks,size(pairIdx,1));

ccTracks = cell(nTracks,1);

marked = false(nTracks);

for iTrack = 1:nTracks
    if ~marked(iTrack)
        ccTracks{iTrack} = rec(iTrack,pairIdx,marked);
    end
end

    function trackInCC = rec(iTrack)
        marked(iTrack) = true;
        
        jdx = find(G(i,:));
        
        for jj = 1:numel(jdx)
            j = jdx(jj);
            
            if (~marked(j))
                
            end
        end
    end

grouping(1:nFrames) = struct('edges',[],'vertices',[]);

for iFrame = 1:nFrames
    % Find which tracks live in iFrame
    inFrameTrackIdx = find(iFirst <= iFrame & iLast >= iFrame);
    
    % Find the pair whose tracks are both in iFrame
    isInFramePair = ismember(pairIdx(:,1), inFrameIdx) & ...
        ismember(pairIdx(:,2), inFrameTrackIdx);
    
    inFramePairIdx = pairIdx(isInFramePair,:);
       
    % Indices of tracks that live in iFrame AND are referred by pairIdx
    trackIdx = sort(unique(inFramePairIdx(:)));
    
    % We need to transform
    %
    % inFramePairIdx = [12 36
    %                   21 44
    %                   26 36]
    %
    % into
    %
    % inFramePairIdx = [1 4
    %                   2 5
    %                   3 4]
    ctr = 1;
    for iTrack = trackIdx
        inFramePairIdx(inFramePairIdx == iTrack) = ctr;
        ctr = ctr + 1;
    end
    
    % the edge of the graph (inFramePairIdx) is ok.
    edges = inFramePairIdx;
    vertices = trackInfos(trackIdx);
    
    % Calculate the connected component and store the CC and the tracks
    % (instead of storing the edges only).
end


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

imagePath = fullfile(movieData.imageDirectory, movieData.channelDirectory{1});
imageFiles = dir([imagePath filesep '*.tif']);
ima = imread(fullfile(imagePath, imageFiles(1).name));

load(fullfile(movieData.particleDetection.directory, ...
    movieData.particleDetection.filename));

X = [featuresInfo(1).xCoord, featuresInfo(1).yCoord];
ind = sub2ind(size(ima),X(:,2),X(:,1));
N = size(ind,1);

[~, T] = steerableFiltering(double(ima),2,2);

Y = [cos(T(ind)), sin(T(ind))];

pair = pcombs(1:N,false);

u0 = X(pair(:,1),:);
u1 = u0 + Y(pair(:,1),:);
v0 = X(pair(:,2),:);
v1 = v0 + Y(pair(:,2),:);

isValid = true(size(pair,1),1);

% euclidian distance [0...+inf]
dE = sqrt(sum((u0 - v0).^2,2));
isValid = isValid & dE <= thE;

% angle between u and v
% dot = abs(sum((u1 - u0) .* (v1 - v0),2));
% dot(dot > 1) = 1;
% dot(dot < -1) = -1;
% dA = acos(dot);
% isValid = isValid & dA <= thA;
% 
% mean distance of u1 and v1 projected on the line (u0,v0)
dP1 = sqrt(1 - (sum((u1-u0) .* (v0-u0),2) ./ dE).^2);
dP2 = sqrt(1 - (sum((v1-v0) .* (v0-u0),2) ./ dE).^2);
dP = dP1 .* dP2;
isValid = isValid & dP <= thP;
% 
% cost = exp(- (dE .* (1/pi) .* dA .* dP));

cost = 1./ dE;

% Build cost matrix
i = pair(isValid,1);
j = pair(isValid,2);
c = cost(isValid);

% Populate the lower triangular part only
D = sparse(j, i, c, N, N, numel(c));

% Compute Maximum Weight Matching
M = maxWeightMatching(D);


% Display result
imshow(ima,[]);
hold on;

B = zeros(N);
ind = sub2ind([N N],M(:,1),M(:,2));
B(ind) = 1;

line(X(:,1),X(:,2),'LineStyle','none', 'Marker', '.', 'Color', 'g');
gplot(B,X,'r');
quiver(X(:,1),X(:,2),Y(:,1),Y(:,2),0,'b');

end