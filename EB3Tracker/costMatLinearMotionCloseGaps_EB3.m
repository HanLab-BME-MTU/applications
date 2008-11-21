function [costMat,nonlinkMarker,indxMerge,numMerge,indxSplit,numSplit,...
    errFlag] = costMatLinearMotionCloseGaps_EB3(trackedFeatInfo,...
    trackedFeatIndx,trackStartTime,trackEndTime,costMatParam,gapCloseParam,...
    kalmanFilterInfo,nnDistLinkedFeat,probDim,movieInfo)
%COSTMATLINEARMOTIONCLOSEGAPS provides a cost matrix for closing gaps and capturing merges/splits using Kalman filter information
%
%SYNOPSIS [costMat,nonlinkMarker,indxMerge,numMerge,indxSplit,numSplit,...
%    errFlag] = costMatLinearMotionCloseGaps(trackedFeatInfo,...
%    trackedFeatIndx,trackStartTime,trackEndTime,costMatParam,gapCloseParam,...
%    kalmanFilterInfo,nnDistLinkedFeat,probDim,movieInfo)
%
%INPUT  trackedFeatInfo: The positions and amplitudes of the tracked
%                        features from linkFeaturesKalman.
%                        Number of rows = number of tracks.
%                        Number of columns = 8*number of frames.
%                        Each row consists of
%                        [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                        in image coordinate system (coordinates in
%                        pixels). NaN is used to indicate time points
%                        where the track does not exist.
%       trackedFeatIndx: Connectivity matrix of features between frames.
%                        Rows indicate continuous tracks, while columns
%                        indicate frames. A track that ends before the
%                        last time point is followed by zeros, and a track
%                        that starts at a time after the first time point
%                        is preceded by zeros.
%       trackStartTime : Starting time of all tracks.
%       trackEndTime   : Ending time of all tracks.
%       costMatParam   : Structure containing variables needed for cost
%                        calculation. Contains the fields:
%             .minSearchRadiusCG:Minimum allowed search radius (in pixels).
%             .maxSearchRadiusCG:Maximum allowed search radius (in pixels).
%                               This value is the maximum search radius
%                               between two consecutive frames as when
%                               linking between consecutive frames. It will
%                               be calcualted for different time gaps
%                               based on the scaling factor of Brownian
%                               motion (expanding it will make use of the
%                               field .timeReachConfB).
%             .brownStdMultCG : Factor multiplying Brownian
%                               displacement std to get search radius.
%                               Vector with number of entries equal to
%                               gapCloseParam.timeWindow (defined below).
%             .linStdMultCG   : Factor multiplying linear motion std to get
%                               search radius. Vector with number of entries
%                               equal to gapCloseParam.timeWindow (defined
%                               below).
%             .timeReachConfB : Time gap for reaching confinement for
%                               2D Brownian motion. For smaller time gaps,
%                               expected displacement increases with
%                               sqrt(time gap). For larger time gaps,
%                               expected displacement increases slowly with
%                               (time gap)^0.1.
%             .timeReachConfL : Time gap for reaching confinement for
%                               linear motion. Time scaling similar to
%                               timeReachConfB above.
%             .lenForClassify : Minimum length of a track to classify it as
%                               directed or Brownian.
%             .maxAngleES     : Maximum allowed angle between two
%                               directions of motion for potential linking
%                               (in degrees).
%             .maxAngleVD     : Max allowed angle between each linear track
%                               and the vector connecting the centers of
%                               the tracks (tests for colinearity)
%             .closestDistScaleCG:Scaling factor of nearest neighbor
%                                 distance.
%             .maxStdMultCG   : Maximum value of factor multiplying
%                               std to get search radius.
%             .ampRatioLimitCG: Minimum and maximum allowed ratio between
%                               the amplitude of a merged feature and the
%                               sum of the amplitude of the two features
%                               making it.
%             .lftCdf         : Lifetime cumulative density function.
%                               Column vector, specifying cdf for
%                               lifetime = 0 to movie length.
%                               Enter [] if cdf is not to be used.
%                               Optional. Default: [].
%             .useLocalDensity: 1 if local density of features is used to expand
%                               their search radius if possible, 0 otherwise.
%             .nnWindow       : Time window to be used in estimating the
%                               nearest neighbor distance of a track at its start
%                               and end.
%       gapCloseParam  : Structure containing variables needed for gap closing.
%                        Contains the fields:
%             .timeWindow : Largest time gap between the end of a track and the
%                           beginning of another that could be connected to it.
%             .tolerance  : Relative change in number of tracks in two
%                           consecutive gap closing steps below which
%                           iteration stops.
%             .mergeSplit : Logical variable with value 1 if the merging
%                           and splitting of trajectories are to be consided;
%                           and 0 if merging and splitting are not allowed.
%       kalmanFilterInfo:Structure array with number of entries equal to
%                        number of frames in movie. Contains the fields:
%             .stateVec   : Kalman filter state vector for each
%                           feature in frame.
%             .stateCov   : Kalman filter state covariance matrix
%                           for each feature in frame.
%             .noiseVar   : Variance of state noise for each
%                           feature in frame.
%             .stateNoise : Estimated state noise for each feature in
%                           frame.
%             .scheme     : 1st column: propagation scheme connecting
%                           feature to previous feature. 2nd column:
%                           propagation scheme connecting feature to
%                           next feature.
%       nnDistLinkedFeat:Matrix indicating the nearest neighbor
%                        distances of features linked together within
%                        tracks.
%       probDim        : Problem dimensionality. 2 (for 2D) or 3 (for 3D).
%       movieInfo      : movieInfo as input to trackCloseGapsKalman. Not
%                        really used in this code, but needed for
%                        compatibility with other cost functions.
%
%OUTPUT costMat       : Cost matrix.
%       nonlinkMarker : Value indicating that a link is not allowed.
%       indxMerge     : Index of tracks that have possibly merged with
%                       tracks that end before the last time points.
%       numMerge      : Number of such tracks.
%       indxSplit     : Index of tracks from which tracks that begin after
%                       the first time point might have split.
%       numSplit      : Number of such tracks.
%       errFlag       : 0 if function executes normally, 1 otherwise.
%
%REMARKS
%
%Khuloud Jaqaman, April 2007

%% Output

costMat = [];
nonlinkMarker = [];
indxMerge = [];
numMerge = [];
indxSplit = [];
numSplit = [];
errFlag = [];

%% Input

%check whether correct number of input arguments was used
if nargin ~= nargin('costMatLinearMotionCloseGaps_EB3')
    disp('--costMatLinearMotionCloseGaps_EB3: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%get cost matrix parameters
minSearchRadius = costMatParam.minSearchRadius;
maxSearchRadius = costMatParam.maxSearchRadius;

brownStdMult = costMatParam.brownStdMult;
timeReachConfB = costMatParam.timeReachConfB;
linStdMult   = costMatParam.linStdMult;
timeReachConfL = costMatParam.timeReachConfL;

lenForClassify = costMatParam.lenForClassify;
useLocalDensity = costMatParam.useLocalDensity;

maxFAngle = (sin(costMatParam.maxFAngle*pi/180))^2;
maxBAngle = (sin(costMatParam.maxBAngle*pi/180))^2;
backVelMultFactor = costMatParam.backVelMultFactor;

nnWindow = costMatParam.nnWindow;
if useLocalDensity
    closestDistScale = 2;
    maxStdMult = 100;
else
    closestDistScale = [];
    maxStdMult = [];
end

if isfield('costMatParam','lftCdf')
    lftCdf = costMatParam.lftCdf;
    oneMinusLftCdf = 1 - lftCdf;
else
    lftCdf = [];
end

%get gap closing parameters
timeWindow = gapCloseParam.timeWindow;
mergeSplit = gapCloseParam.mergeSplit;

%make sure that timeReachConfB and timeReachConfL are <= timeWindow
timeReachConfB = min(timeReachConfB,timeWindow);
timeReachConfL = min(timeReachConfL,timeWindow);

%find the number of tracks to be linked and the number of frames in the movie
[numTracks,numFrames] = size(trackedFeatInfo);
numFrames = numFrames / 8;

%list the tracks that start and end in each frame
tracksPerFrame = repmat(struct('starts',[],'ends',[]),numFrames,1);
for iFrame = 1 : numFrames
    tracksPerFrame(iFrame).starts = find(trackStartTime == iFrame); %starts
    tracksPerFrame(iFrame).ends = find(trackEndTime == iFrame); %ends
end

%% Gap closing

%get the x,y-coordinates and amplitudes at the starts of tracks
coordStart = zeros(numTracks,probDim);
ampStart   = zeros(numTracks,1);
for iTrack = 1 : numTracks
    coordStart(iTrack,:) = trackedFeatInfo(iTrack,...
        (trackStartTime(iTrack)-1)*8+1:(trackStartTime(iTrack)-1)*8+probDim);
    ampStart(iTrack) = trackedFeatInfo(iTrack,(trackStartTime(iTrack)-1)*8+4);
end

%get the x,y-coordinates and amplitudes at the ends of tracks
coordEnd = zeros(numTracks,probDim);
ampEnd   = zeros(numTracks,1);
for iTrack = 1 : numTracks
    coordEnd(iTrack,:) = trackedFeatInfo(iTrack,...
        (trackEndTime(iTrack)-1)*8+1:(trackEndTime(iTrack)-1)*8+probDim);
    ampEnd(iTrack) = trackedFeatInfo(iTrack,(trackEndTime(iTrack)-1)*8+4);
end

%determine the types, velocities, noise stds and centers of all tracks
[trackType,xyzVel,noiseStd,trackCenter] = estimTrackTypeParamLM(...
    trackedFeatIndx,trackedFeatInfo,kalmanFilterInfo,lenForClassify,probDim);

percentUnclass = sum(isnan(trackType))/length(trackType);
percentBrown = sum((trackType==0))/length(trackType);
percentLinear = sum((trackType==1))/length(trackType);

save('trackTypes', 'trackType','percentUnclass','percentBrown','percentLinear')

% get vel from kalman filter for tracks that were brownian or unclassified
for iTrack=1:numTracks
    if trackType(iTrack)~=1
        xyzVel(iTrack,:) = kalmanFilterInfo(trackEndTime(...
            iTrack)).stateVec(trackedFeatIndx(iTrack,...
            trackEndTime(iTrack)),probDim+1:2*probDim);
    end
end


trackTypeTemp = ones(size(trackType)); % make all tracks artificially LINEAR

%find the 80th percentile of the noise standard deviation in order to use
%that for undetermined tracks (after removing std = 1 which indicates the
%simple initialization conditions
noiseStdAll = noiseStd(noiseStd ~= 1);
undetBrownStd = prctile(noiseStdAll,80);

%determine the search areas of all tracks - here we trick the function by
%considering all tracks as linear. that is, we have gotten xyzVel and
%reassigned trackType to 1 for all.  previously, this function was used to
%get a search rectangle.  now we just use longVecS and longVecE, which
%determine the distance criterion for whether a start falls within the
%end's cone...
[longVecSAll,longVecEAll,shortVecSAll,shortVecEAll,shortVecS3DAll,...
    shortVecE3DAll] = getAveDispEllipseAll(xyzVel,noiseStd,trackTypeTemp,...
    undetBrownStd,timeWindow,brownStdMult,linStdMult,timeReachConfB,...
    timeReachConfL,minSearchRadius,maxSearchRadius,useLocalDensity,...
    closestDistScale,maxStdMult,nnDistLinkedFeat,nnWindow,...
    trackStartTime,trackEndTime,probDim);

%find all pairs of ends and starts that can potentially be linked
%determine this by looking at time gaps between ends and starts
%and by looking at the distance between pairs
indxEnd2 = [];
indxStart2 = [];

%get the absolute upper limit of acceptable displacements in one frame
%as the maximum of (maximum velocity multiplied by probDim*linStdMult(1),
%maxSearchRadiusCG)
maxDispAllowed = max(max(xyzVel(:)) * probDim * linStdMult(1), ...
    maxSearchRadius);

%go over all frames until the one before last
for iFrame = 1 : numFrames - 1

    %find tracks that end in this frame
    endsToConsider = tracksPerFrame(iFrame).ends;

    for jFrame = iFrame + 1 : min(iFrame+timeWindow,numFrames)

        %find tracks that start in this frame
        startsToConsider = tracksPerFrame(jFrame).starts;

        %calculate the distance between ends and starts
        dispMat2 = createDistanceMatrix(coordEnd(endsToConsider,:),...
            coordStart(startsToConsider,:));

        %find possible pairs
        [indxEnd3,indxStart3] = find(dispMat2 <= maxDispAllowed * sqrt(jFrame-iFrame));
        if size(indxEnd3,1) == 1
            indxEnd3 = indxEnd3';
            indxStart3 = indxStart3';
        end

        %add them to the list of possible pairs
        indxEnd2 = [indxEnd2; endsToConsider(indxEnd3)];
        indxStart2 = [indxStart2; startsToConsider(indxStart3)];

    end %(for jFrame = iFrame + 1 : iFrame + timeWindow)

end %(for iFrame = 1 : numFrames)

%get total number of pairs
numPairs = length(indxEnd2);

%clear variables from memory
clear dispMat2 maxDispAllowed

%reserve memory for cost matrix vectors
indx1 = zeros(numPairs,1); %row number in cost matrix
indx2 = zeros(numPairs,1); %column number in cost matrix

%for cost
pairDist = zeros(numPairs,1);
pairAngle1 = zeros(numPairs,1);
pairAngle2 = zeros(numPairs,1);
pairTrackType = zeros(numPairs,2); %fyi, track types

%put time scaling of linear motion in a vector
timeScalingLin = ones(timeWindow,1);

%put time scaling of Brownian motion in a vector
timeScalingBrown = ones(timeWindow,1);

%go over all possible pairs of starts and ends
for iPair = 1 : numPairs

    %get indices of start and end and the time gap between them
    iStart = indxStart2(iPair);
    iEnd = indxEnd2(iPair);

    %determine the time gap between them
    timeGap = trackStartTime(iStart) - trackEndTime(iEnd);

    %get the types of the two tracks
    trackTypeS = trackType(iStart);
    trackTypeE = trackType(iEnd);

    %calculate vectors whose lengths correspond to supposed distance
    %traveled within the time gap, in either forward (F) or backward (B) direction
    longVecSF = longVecSAll(:,timeGap,iStart);
    longVecSB = -backVelMultFactor.*longVecSF;
    longVecEF = longVecEAll(:,timeGap,iEnd);
    longVecEB = -backVelMultFactor.*longVecEF;

    %calculate the magnitudes of the long vectors of both start and end
    longVecMagSF = norm(longVecSF);
    longVecMagSB = norm(longVecSB);
    longVecMagEF = norm(longVecEF);
    longVecMagEB = norm(longVecEB);

    %calculate the vector connecting the end of track iEnd to the
    %start of track iStart and compute its magnitude
    dispVec = coordStart(iStart,:) - coordEnd(iEnd,:); % kat (flipped direction)
    dispVecMag = norm(dispVec);

    
    
    %cosAngle between the end and start velocity vectors
    %if positive, the tracks point in same general direction
    %if negative, don't allow a link because comets on the same MT must
    %travel in the same general direction
    cosAngleEF_SF = (longVecEF' * longVecSF) / (longVecMagEF * longVecMagSF);
    
    %cosAngle between the end forward vector and the displacement vector
    cosAngleEF_D = (longVecEF' * dispVec') / (longVecMagEF * dispVecMag);

    %cosAngle between the end backward vector and the displacement vector
    cosAngleEB_D = (longVecEB' * dispVec') / (longVecMagEB * dispVecMag);

    
    if cosAngleEF_SF < 0  
        %tracks point in opposite directions
        possibleLink = 0;
    else
        %tracks point in same general direction
        if cosAngleEF_D > 0 
            %check whether start meets angle and distance criteria for
            %end's forward cone
            possibleLink = cosAngleEF_D >= cos(maxFAngle) && dispVecMag <= longVecMagEF;         
        else
            %check whether start meets angle and distance criteria for
            %end's backward cone
            possibleLink = cosAngleEB_D >= cos(maxBAngle) && dispVecMag <= longVecMagEB; 
        end
    end

    %calculate square sine of angle between each motion direction vector
    %and the displacement vector...both are close to 0 when vectors are
    %parallel
    sin2AngleEF_D = 1 - (dispVec * longVecEF / (dispVecMag * longVecMagEF))^2;
    sin2AngleSF_D = 1 - (dispVec * longVecSF / (dispVecMag * longVecMagSF))^2;

    if possibleLink
        %get costs from distance and angles
        pairDist(iPair) = dispVecMag;
        pairAngle1(iPair) = sin2AngleEF_D;
        pairAngle2(iPair) = sin2AngleSF_D;
        pairTrackType(iPair,:) = [trackTypeE trackTypeS];

        %specify the location of this pair in the cost matrix
        indx1(iPair) = iEnd;   %row number
        indx2(iPair) = iStart; %column number
    end

end %(for iPair = 1 : numPairs)

%keep only pairs that turned out to be possible
possiblePairs = find(indx1 ~= 0);
indx1 = indx1(possiblePairs);
indx2 = indx2(possiblePairs);
pairDist = pairDist(possiblePairs);
pairAngle1 = pairAngle1(possiblePairs);
pairAngle2 = pairAngle2(possiblePairs);
pairTrackType = pairTrackType(possiblePairs,:);

% calculate the costs
costDist = (pairDist - mean(pairDist)).^2/(std(pairDist)).^2;
costAngle1 = (pairAngle1 - mean(pairAngle1)).^2/(std(pairAngle1)).^2;
costAngle2 = (pairAngle2 - mean(pairAngle2)).^2/(std(pairAngle2)).^2;

%cost = costAngle1 + costAngle2;
cost = pairAngle1 + pairAngle2; % use only angle info - will optimize straightest link
clear possiblePairs

%% Merging and splitting

%define some merging and splitting variables
numMerge  =  0; %index counting merging events
indxMerge = []; %vector storing merging track number
altCostMerge = []; %vector storing alternative costs of not merging
numSplit  =  0; %index counting splitting events
indxSplit = []; %vector storing splitting track number
altCostSplit = []; %vector storing alternative costs of not splitting
%create cost matrix without births and deaths
numEndSplit = numTracks;
numStartMerge = numTracks;
costMat = sparse(indx1,indx2,cost,numEndSplit,numStartMerge);

%% Append cost matrix to allow births and deaths ...

%determine the cost of birth and death
costBD = prctile(cost,90);

%get the cost for the lower right block
costLR = min(min(min(costMat))-1,-1);

% % %create cost matrix that allows for births and deaths
costMat = [costMat ... %costs for links (gap closing + merge/split)
    spdiags([costBD*ones(numTracks,1); altCostSplit],0,numEndSplit,numEndSplit); ... %costs for death
    spdiags([costBD*ones(numTracks,1); altCostMerge],0,numStartMerge,numStartMerge) ...  %costs for birth
    sparse(indx2,indx1,costLR*ones(length(indx1),1),numStartMerge,numEndSplit)]; %dummy costs to complete the cost matrix

%determine the nonlinkMarker
nonlinkMarker = min(floor(full(min(min(costMat))))-5,-5);


%% ~~~ the end ~~~
