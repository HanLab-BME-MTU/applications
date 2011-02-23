function movieData = getMoviePairTracks(movieData,minOverlap,timeGap,maxEuclidianDist,kSigma,nLevels,alpha,probBinSize,batchMode)

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
% kSigma:               cutoff in number of standard deviations
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

imagePath = fullfile(movieData.imageDirectory, movieData.channelDirectory{1});
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
% x11 dx11 y11 dy11 A11 dA11 Sx11 dSx11 Sy11 dSy11 T11 dT11 C11 dC11 varE11
% x12 dx12 ... (2nd point of the 1st track)
% ...
% x1n dx1n ... (nth point of the 1st track)
% x21 dx21 ... (1st point of the 2nd track)
% ...
params = nan(numel(tracksFeatIndx),15);

for iFrame = 1:nFrames
    ind = frameTracks == iFrame & tracksFeatIndx ~= 0;
    xCoord = [featuresInfo(iFrame).xCoord];
    yCoord = [featuresInfo(iFrame).yCoord];
    amp = [featuresInfo(iFrame).amp];
    sX = [featuresInfo(iFrame).stdAlong];
    sY = [featuresInfo(iFrame).stdAside];
    theta = [featuresInfo(iFrame).theta];
    bkg = [featuresInfo(iFrame).bkg];
    varError = [featuresInfo(iFrame).varError];
    
    params(ind,1:2) = xCoord(tracksFeatIndx(ind),:);
    params(ind,3:4) = yCoord(tracksFeatIndx(ind),:);
    params(ind,5:6) = amp(tracksFeatIndx(ind),:);
    params(ind,7:8) = sX(tracksFeatIndx(ind),:);
    params(ind,9:10) = sY(tracksFeatIndx(ind),:);
    params(ind,11:12) = theta(tracksFeatIndx(ind),:);
    params(ind,13:14) = bkg(tracksFeatIndx(ind),:);
    params(ind,15) = varError(tracksFeatIndx(ind));
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
    params(gacombIdx,13:14) = interp1(borderIdx, params(borderIdx,13:14), gacombIdx);
    params(gacombIdx,15) = interp1(borderIdx, params(borderIdx,15), gacombIdx);
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

% % (DEBUG) Save pair tracks per frame
% trackLabels = computeLabel(validTrackPairIdx, nTracks); %#ok<NASGU>
% save(fullfile(movieData.pairTracks.directory, ...
%     'classifiedTracksPerDistance.mat'),'tracksFinal', 'trackLabels');

% 6) Hierarchical Track Clustering

% CC: a cell array of connected components. Each connected component
% contains a list of track indices. At the beginning, Each connected
% component contains one track only.
CC = arrayfun(@(x) {x}, 1:nTracks);
nCC = numel(CC);

% CCmodels: is a struct array containing the following fields:
% .iFirst:     array of nCCx1 corresponding to the birth frame of each CC
% .iLast:      array of nCCx1 corresponding to the death frame of each CC
% .type:       array of nCCx1 taking the following values:
%                 1 == anisotropic model
%                 2 == segment model
% .paramsInd:  array of nCCx1 containing indices is .modelParams{*}
%              array. The index in the cell array .modelParams is
%              given by .modelType.
% .params:     cell array of size 2 containing model parameter of
%              anisotropic model (index 1) and segment models (index 2)
%                if modelType == 1:
%                 x1 dx1 y1 dy1 A1 dA1 Sx1 dSx1 Sy1 dSy1 T1 dT1 C1 dC1 var1
%                 x2 dx2 y2 dy2 A2 dA2 Sx2 dSx2 Sy2 dSy2 T2 dT2 C2 dC2 var2
%                 ...
%                 xN dxN yN ... varN
%
%                if modelType == 2:
%                 x1 dx1 y1 dy1 A1 dA1 l1 dl1 s1 ds1 T1 dT1 C1 dC1 var1
%                 x2 dx2 y2 dy2 A2 dA2 l2 dl2 s2 ds2 T2 dT2 C2 dC2 var2
%                 ...
%                 xN dxN yN ... varN

CCmodels.iFirst = iFirst;
CCmodels.iLast = iLast;
CCmodels.type = ones(1, nCC);
CCmodels.paramsFirst = first;
CCmodels.paramsLast = last;
CCmodels.params = cell(2,1);
CCmodels.params{1} = params;

% CCpairIdx: pair of CC indices for which we need a weight. We initialize
% CCpairIdx with validTrackPairIdx
CCpairIdx = validTrackPairIdx;
nCCPairs = size(CCpairIdx,1);

initFunc = {...
    @initPairParams00,...
    @initPairParams01,...
    @(params1,params2) initPairParams01(params2,params1),...
    @initPairParams11};
    
for iLevel = 1:nLevels

    % Compute model parameters of each pair of CC
    iFirst1 = CCmodels.iFirst(CCpairIdx(:,1));
    iFirst2 = CCmodels.iFirst(CCpairIdx(:,2));
    iLast1 = CCmodels.iLast(CCpairIdx(:,1));
    iLast2 = CCmodels.iLast(CCpairIdx(:,2));
    
    iFirstPair = max(iFirst1,iFirst2);
    iLastPair = min(iLast1,iLast2);
    overlap = iLastPair - iFirstPair + 1;
    offset1 = iFirstPair - iFirst1;
    offset2 = iFirstPair - iFirst2;
        
    outputFirst = cumsum([1 overlap(1:end-1)]);
    
    % Compute model parameters of all pairs.
    % 8 = 7 parameters (x,y,A,l,sigma,theta,C) + varError
    CCpairParams = zeros(sum(overlap), 8);
    
    % iFunc == 1: modelType1 == 1, modelType2 == 1
    % iFunc == 2: modelType1 == 1, modelType2 == 2
    % iFunc == 3: modelType1 == 2, modelType2 == 1
    % iFunc == 4: modelType1 == 2, modelType2 == 2
    for iFunc = 1:4
        modelType1 = ceil(iFunc/2);
        modelType2 = mod(iFunc-1,2)+1;
        
        isFunc = ...
            CCmodels.type(CCpairIdx(:,1)) == modelType1 & ...
            CCmodels.type(CCpairIdx(:,2)) == modelType2;
        
        inputInd1 = cell2mat(arrayfun(@(x,y) x:(x+y-1), ...
            CCmodels.paramsFirst(CCpairIdx(isFunc,1)) + offset1(isFunc), ...
            overlap(isFunc), 'UniformOutput',false));
        
        inputInd2 = cell2mat(arrayfun(@(x,y) x:(x+y-1), ...
            CCmodels.paramsFirst(CCpairIdx(isFunc,2)) + offset2(isFunc), ...
            overlap(isFunc), 'UniformOutput', false));
        
        outputInd = cell2mat(arrayfun(@(x,y) x:(x+y-1), ...
            outputFirst(isFunc), overlap(isFunc), ...
            'UniformOutput', false));
        
        CCpairParams(outputInd,1:end-1) = initFunc{iFunc}(...
            CCmodels.params{modelType1}(inputInd1,:), ...
            CCmodels.params{modelType2}(inputInd2,:));
    end
    
    % DEBUG: save CCpairParams as featureInfo
    %segments = cell(nFrames,1);
    
    if ~batchMode
        h = waitbar(0,['Please wait, computing pair models at level' ...
            num2str(iLevel) '...']);
    end

    % Compute pair model
    for iFrame = 1:nFrames
        % Read image
        ima = double(imread(fullfile(imagePath, imageFiles(iFrame).name)));
        [ny nx] = size(ima);
        
        % Find which pair of CC is living at frame iFrame
        isInFrame = iFirstPair <= iFrame & iLastPair >= iFrame;

        ind = outputFirst(isInFrame) + iFrame - iFirstPair(isInFrame);
        
        %segments{iFrame} = CCpairParams(ind, :);
        
        initParams = CCpairParams(ind,1:end-1);
        nInitParams = size(initParams);
        
        for iParam = 1:nInitParams
            
            cellParams = num2cell(initParams(iParam,:),1);
            [x0 y0 A l sigma theta C] = cellParams{:};
            [xRange,yRange,nzIdx] = ...
                segment2DSupport(x0,y0,l,sigma,theta,[ny nx],kSigma);

            crop = ima(yRange,xRange);
            mask = false(size(crop));
            mask(nzIdx) = true;
            crop(~mask) = NaN;            
            
            dx = x0 - xRange(floor((numel(xRange)+1)/2));
            dy = y0 - yRange(floor((numel(yRange)+1)/2));
            
            [prm, ~, ~, R] = fitSegment2D(crop, [dx, dy, A, l, sigma, theta, C], 'AlC');
            
            CCpairParams(ind(iParam),1:end-1) = prm;
            CCpairParams(ind(iParam),8) = (1/(numel(nzIdx)-1)) * sum(R(nzIdx).^2);
        end
        
        if ~batchMode && ishandle(h)
            waitbar(iFrame/nFrames, h)
        end
    end
    
    if ~batchMode && ishandle(h)
        close(h);
    end

    %save(fullfile(movieData.pairTracks.directory, 'CCpairParams.mat'), 'segments');

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
params1 = num2cell(params1,1);
params2 = num2cell(params2,1);
[x1,y1,A1,sx1,sy1,~,C1,~] = params1{1:2:end};
[x2,y2,A2,sx2,sy2,~,C2,~] = params2{1:2:end};

t0 = atan2(y2-y1,x2-x1);
ct = cos(t0);
st = sin(t0);
ex1 = x1 - sx1 .* ct;
ey1 = y1 - sx1 .* st;
ex2 = x2 + sx2 .* ct;
ey2 = y2 + sx2 .* st;
x0 = .5 * (ex1+ex2);
y0 = .5 * (ey1+ey2);
A0 = .5 * (A1+A2);
l0 = sqrt((ex2-ex1).^2 + (ey2-ey1).^2);
s0 = .5 * (sy1+sy2);
C0 = .5 * (C1+C2);
pairParams = [x0 y0 A0 l0 s0 t0 C0];

function pairParams = initPairParams01(params1,params2)
% TODO !!!
pairParams = zeros(size(params1,1),7);

function pairParams = initPairParams11(params1,params2)
% TODO !!!
pairParams = zeros(size(params1,1),7);

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

