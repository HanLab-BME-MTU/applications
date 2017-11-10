function  getMoviePairTracks(featuresInfo, trackFinal, sigmaPSF, kSigma, nFrames,imSize,pixelSize, distTransPath, outputPath, varargin)

% minLifetime:     is the minimal number of frames every tracks must last.
%
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

% BEGIN
%movieData.pairTracks.status = 0;

% Parse input parameters
%checkMovieData = @(movieData) ...
%    checkMovieDistanceTransform(movieData) && ...
%    checkMovieParticleDetection(movieData) && ...
%    checkMovieParticleTracking(movieData);
check = (@(x) isstruct(x) | ischar(x)); 
ip = inputParser;
ip.CaseSensitive = false;

ip.addRequired('featuresInfo', check);
ip.addRequired('trackFinal',   check);
ip.addRequired('sigmaPSF',     @isscalar);
ip.addRequired('kSigma',       @isscalar);
ip.addRequired('nFrames',      @isscalar);
ip.addRequired('imSize',       @isvector);
ip.addRequired('pixelSize',    @isscalar);
ip.addRequired('distTransPath',@ischar);
ip.addRequired('outputPath',   @ischar);

ip.addOptional('minLifetime', 1,    @isscalar);
ip.addOptional('maxDistance', 2000, @isscalar);
ip.addOptional('minOverlap',  1,    @isscalar);
ip.addOptional('bandWidth',   1000, @isscalar);
ip.addOptional('minDistance', 350,  @isscalar);
ip.addOptional('alpha',       0.05, @isscalar);

if ischar(featuresInfo)
    aux = load(featuresInfo);
    clear featuresInfo
    featuresInfo = aux.featuresInfo;
end

if ischar(trackFinal)
    aux = load(trackFinal);
    clear trackFinal
    trackFinal = aux.trackFinal;
end

ip.parse(featuresInfo, trackFinal, sigmaPSF, kSigma, nFrames,imSize,pixelSize, distTransPath, outputPath, varargin{:});

minLifetime = ip.Results.minLifetime;
maxDistance = ip.Results.maxDistance;
minOverlap  = ip.Results.minOverlap;
bandWidth   = ip.Results.bandWidth;
minDistance = ip.Results.minDistance;
alpha       = ip.Results.alpha;


% Get all track parameters
[tracksFinal, allFeatures, tFirst, lifetime] = getAllFeatures(featuresInfo,trackFinal, minLifetime); %#ok<ASGLU>
%tracksFinal is updated in this function
%


nTracks = numel(tFirst);
tLast = tFirst + lifetime - 1;

% pFirst and pLast are indexing allFeatures
pLast = cumsum(lifetime);
pFirst = pLast-lifetime+1;

% Get initial pair track candidates
E = getInitialPairTracks(nFrames,imSize,pixelSize, distTransPath, allFeatures, tFirst, lifetime, ...
    maxDistance, minOverlap, bandWidth, minDistance);

% Define the set of connected components
% numel(tFirst) == nTracks
CC = arrayfun(@(t) {t}, (1:nTracks)');
nCC = numel(CC);

% Start iteration
iter = 0;

while size(E,1)
    iter = iter + 1;
    
    fprintf(1, 'Iterative Track Clustering level %d\n', iter);
    fprintf(1, '\tNumber of connected components:\t%d\n', nCC);
    fprintf(1, '\tNumber of pair candidates:\t%d\n', size(E,1));   
    
    % Compute the first and last frame of each CC
    tFirstCC = cellfun(@(trackIdx) min(tFirst(trackIdx)), CC); % first frame of CC
    tLastCC  = cellfun(@(trackIdx) max(tLast(trackIdx)), CC);   % last frame of CC
    
    % Compute the overlap between CC
    tOverlapFirst = max(tFirstCC(E(:,1)), tFirstCC(E(:,2)));
    tOverlapFirst = arrayfun(@(x) {x}, tOverlapFirst);
    
    tOverlapLast = min(tLastCC(E(:,1)), tLastCC(E(:,2)));
    tOverlapLast = arrayfun(@(x) {x}, tOverlapLast);

    % For each track in the first/second CC, compute the first feature
    % index that appears at tOverlapFirst
    pFirstCC1 = cellfun(@(trackIdx, tOverlapFirst) ...
        pFirst(trackIdx) + tOverlapFirst - tFirst(trackIdx), ...
        CC(E(:,1)), tOverlapFirst, 'UniformOutput', false);
    pFirstCC2 = cellfun(@(trackIdx, tOverlapFirst) ...
        pFirst(trackIdx) + tOverlapFirst - tFirst(trackIdx), ...
        CC(E(:,2)), tOverlapFirst, 'UniformOutput', false);
    
    % For each track in the first/second CC, compute its lifetime between
    % tOverlapFirst and tOverlapLast
    trackLifetimeInOverlapCC1 = cellfun(@(trackIdx, tOverlapFirst, tOverlapLast) ...
        min(tLast(trackIdx), tOverlapLast) - max(tFirst(trackIdx), ...
        tOverlapFirst) + 1, CC(E(:,1)), tOverlapFirst, tOverlapLast, ...
        'UniformOutput', false);
    trackLifetimeInOverlapCC2 = cellfun(@(trackIdx, tOverlapFirst, tOverlapLast) ...
        min(tLast(trackIdx), tOverlapLast) - max(tFirst(trackIdx), ...
        tOverlapFirst) + 1, CC(E(:,2)), tOverlapFirst, tOverlapLast, ...
        'UniformOutput', false);
    
    % Gather every feature of each track in the first/second CC between
    % tOverlapFirst and tOverlapLast
    allFeaturesCC1 = cellfun(@(aa,bb) arrayfun(@(a,b) ...
        allFeatures(a:a+b-1, [1 2 4 6]), aa, bb, 'UniformOutput', false), ...
        pFirstCC1, trackLifetimeInOverlapCC1, 'UniformOutput', false);
    
    allFeaturesCC1 = cellfun(@(c) vertcat(c{:}), allFeaturesCC1, ...
        'UniformOutput', false);
    
    allFeaturesCC2 = cellfun(@(aa,bb) arrayfun(@(a,b) ...
        allFeatures(a:a+b-1, [1 2 4 6]), aa, bb, 'UniformOutput', false), ...
        pFirstCC2, trackLifetimeInOverlapCC2, 'UniformOutput', false);
    
    allFeaturesCC2 = cellfun(@(c) vertcat(c{:}), allFeaturesCC2, ...
        'UniformOutput', false);

    % Compute the number of features per each first/second CC
    nFeatures1 = cellfun(@(a) size(a,1), allFeaturesCC1);    
    nFeatures2 = cellfun(@(a) size(a,1), allFeaturesCC2);
   
    % Compute model of the each first/second CC
    [~, res1] = getSegmentModels(allFeaturesCC1);
    [~, res2] = getSegmentModels(allFeaturesCC2);
    
    % Compute model of each pair of CC
    allPairFeatures = cellfun(@(c1,c2) [c1;c2], allFeaturesCC1, ...
        allFeaturesCC2, 'UniformOutput', false);
    
    [~, resPair] = getSegmentModels(allPairFeatures);    
    
    % Compute BIC of the split model versus merged model
    N = 2 * (nFeatures1 + nFeatures2);
    varSplit = cellfun(@(res1,res2) var([res1; res2], 1), res1, res2);
    varMerge = cellfun(@(res) var(res, 1), resPair);

    % If a variance is too small due to too small number of points, we put
    % a lower band of .25. This should be avoided by introducing in the fit
    % the localization error of each point (see regression note) so that
    % even if there are only 2 points to make the regression, the residual
    % won't be 0. The rational behind .25 is that the error of localization
    % is about half a pixel.
    bicSplit = N .* log(max(varSplit,.25)) + 8 * log(N);
    bicMerge = N .* log(max(varMerge,.25)) + 4 * log(N);

    % Define score of each pair
    W = bicSplit - bicMerge;

    % Define which pair is valid
    %isGoodOfFit = cellfun(@(res) ~kstest(res ./ std(res, 1), [], alpha), ...
    %    resPair);
    
    % disable invalid pairs
    ind = find(W >= 0 & varMerge <= (kSigma * sigmaPSF)^2);
    
    if isempty(ind)
        break;
    end
    
    M = maxWeightedMatching(nCC, E(ind,:), W(ind));
        
    % Merge E, CC and update nCC
    EM                  = E(ind(M),:);
    ccInd               = 1:nCC;
    ccInd(EM(:,2))      = EM(:,1);
    ccIndUnique         = unique(ccInd);
    values              = 1:max(ccIndUnique);
    values(ccIndUnique) = 1:numel(ccIndUnique);
    ccInd               = values(ccInd);

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
end

% Clean up CC
[CC allFeatures] = cleanUpCC(nFrames,imSize,pixelSize,distTransPath, CC, allFeatures, tFirst, tLast, ...
    pFirst, bandWidth, alpha);
nCC = numel(CC);

% Save the labeled tracks
trackLabels = zeros(nTracks,1);
for iCC = 1:nCC
    trackLabels(CC{iCC}) = iCC;
end

filename = fullfile(outputPath, 'ClassifiedTracks.mat');
save(filename, 'tracksFinal', 'trackLabels');

% Save the segment models
filename = fullfile(outputPath, 'ClassifiedSegments.mat');
saveCCTracks(nFrames, CC, allFeatures, tFirst, tLast, pFirst, filename);