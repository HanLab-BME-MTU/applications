function [costMat,nonlinkMarker,indxMerge,numMerge,indxSplit,numSplit,...
    errFlag] = makiTrackCostMatCloseGaps(trackedFeatInfo,...
    trackedFeatIndx,trackStartTime,trackEndTime,costMatParam,gapCloseParam,...
    kalmanFilterInfo,nnDistLinkedFeat,probDim,movieInfo)
%MAKITRACKCOSTMATCLOSEGAPS provides a cost matrix for closing gaps in tracking HeLa kinetochores
%
%SYNOPSIS [costMat,nonlinkMarker,indxMerge,numMerge,indxSplit,numSplit,...
%    errFlag] = makiTrackCostMatCloseGaps(trackedFeatInfo,...
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
%                               3 entries for inlier, unaligned and
%                               lagging kinetochores.
%             .maxSearchRadiusCG:Maximum allowed search radius (in pixels).
%                               3 entries for inlier, unaligned and
%                               lagging kinetochores.
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
%             .maxAngle       : Maximum allowed angle between two
%                               directions of motion for potential linking
%                               (in degrees).
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
%                        compatibility with other cost functions. This
%                        code needs the kinType field of movieInfo.
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
%Khuloud Jaqaman, July 2008

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
if nargin ~= nargin('makiTrackCostMatCloseGaps')
    disp('--makiTrackCostMatCloseGaps: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%get cost matrix parameters
linearMotion = costMatParam.linearMotion;
minSearchRadius = costMatParam.minSearchRadius;
maxSearchRadius = costMatParam.maxSearchRadius;
brownStdMult = costMatParam.brownStdMult;
timeReachConfB = costMatParam.timeReachConfB;
lenForClassify = costMatParam.lenForClassify;
useLocalDensity = costMatParam.useLocalDensity;
linStdMult   = costMatParam.linStdMult;
timeReachConfL = costMatParam.timeReachConfL;
sin2AngleMax = (sin(costMatParam.maxAngleVV*pi/180))^2;
sin2AngleMaxVD = 1;
nnWindow = costMatParam.nnWindow;
if useLocalDensity
    closestDistScale = 2;
    maxStdMult = 100;
else
    closestDistScale = [];
    maxStdMult = [];
end
minAmpRatio = costMatParam.ampRatioLimit(1);
maxAmpRatio = costMatParam.ampRatioLimit(2);

if isfield('costMatParam','lftCdf')
    lftCdf = costMatParam.lftCdf;
    oneMinusLftCdf = 1 - lftCdf;
else
    lftCdf = [];
end

if isfield(costMatParam,'gapPenalty') && ~isempty(costMatParam.gapPenalty)
    gapPenalty = costMatParam.gapPenalty;
else
    gapPenalty = 1;
end

%get gap closing parameters
timeWindow = gapCloseParam.timeWindow;
mergeSplit = gapCloseParam.mergeSplit;

%make sure that timeReachConfB and timeReachConfL are <= timeWindow
timeReachConfB = min(timeReachConfB,timeWindow);
timeReachConfL = min(timeReachConfL,timeWindow);

%find the number of track segments to be linked and the number of frames in the movie
[numTracks,numFrames] = size(trackedFeatInfo);
numFrames = numFrames / 8;

%list the track segments that start and end in each frame
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

%determine the kinetochore type at the start and end of each track segment
trackStartKinType = NaN(numTracks,1);
trackEndKinType = trackStartKinType;
for iTrack = 1 : numTracks
    trackStartKinType(iTrack) = movieInfo(trackStartTime(iTrack)).kinType(...
        trackedFeatIndx(iTrack,trackStartTime(iTrack)));
    trackEndKinType(iTrack) = movieInfo(trackEndTime(iTrack)).kinType(...
        trackedFeatIndx(iTrack,trackEndTime(iTrack)));
end

%determine the types, velocities, noise stds and centers of all tracks
[trackType,xyzVel,noiseStd,trackCenter] = estimTrackTypeParamLM(...
    trackedFeatIndx,trackedFeatInfo,kalmanFilterInfo,lenForClassify,probDim);

%if by chance some tracks are labeled linear when linearMotion=0, make them
%Brownian
if linearMotion ~= 1
    trackType(trackType==1) = 0;
end

%find the 80th percentile of the noise standard deviation in order to use
%that for undetermined tracks (after removing std = 1 which indicates the
%simple initialization conditions
noiseStdAll = noiseStd(noiseStd ~= 1);
undetBrownStd = prctile(noiseStdAll,80);

%determine the search areas of all tracks
[longVecSAll,longVecEAll,shortVecSAll,shortVecEAll,shortVecS3DAll,...
    shortVecE3DAll] = makiTrackGetAveDispEllipseAll(xyzVel,noiseStd,trackType,...
    undetBrownStd,timeWindow,brownStdMult,linStdMult,timeReachConfB,...
    timeReachConfL,minSearchRadius,maxSearchRadius,useLocalDensity,...
    closestDistScale,maxStdMult,nnDistLinkedFeat,nnWindow,...
    trackStartTime,trackEndTime,probDim,trackStartKinType,trackEndKinType);

%find all pairs of ends and starts that can potentially be linked
%determine this by looking at time gaps between ends and starts
%and by looking at the distance between pairs
indxEnd2 = [];
indxStart2 = [];

%get the absolute upper limit of acceptable displacements in one frame
%as the maximum of (maximum velocity multiplied by probDim*linStdMult(1),
%maximum maxSearchRadiusCG)
maxDispAllowed = max(max(xyzVel(:)) * probDim * linStdMult(1), ...
    max(maxSearchRadius));

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
cost  = zeros(numPairs,1); %cost value

%put time scaling of linear motion in a vector
% timeScalingLin = [(1:timeReachConfL) timeReachConfL * ...
%     (2:timeWindow-timeReachConfL+1).^0.02];
timeScalingLin = ones(timeWindow,1);

%put time scaling of Brownian motion in a vector
% timeScalingBrown = [(1:timeReachConfB) timeReachConfB * ...
%     (2:timeWindow-timeReachConfB+1).^0.02];
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

    %determine the search area of track iStart
    longVecS = longVecSAll(:,timeGap,iStart);
    shortVecS = shortVecSAll(:,timeGap,iStart);

    %determine the search area of track iEnd
    longVecE = longVecEAll(:,timeGap,iEnd);
    shortVecE = shortVecEAll(:,timeGap,iEnd);

    %calculate the magnitudes of the long and short search vectors
    %of both start and end
    longVecMagS = norm(longVecS);
    shortVecMagS = norm(shortVecS);
    longVecMagE = norm(longVecE);
    shortVecMagE = norm(shortVecE);

    %calculate the vector connecting the end of track iEnd to the
    %start of track iStart and compute its magnitude
    dispVec = coordEnd(iEnd,:) - coordStart(iStart,:);
    dispVecMag = norm(dispVec);

    %project the connecting vector onto the long and short vectors
    %of track iStart and take absolute value
    projStartLong = abs(dispVec * longVecS) / longVecMagS;
    projStartShort = abs(dispVec * shortVecS) / shortVecMagS;
    
    %project the connecting vector onto the long and short vectors
    %of track iEnd and take absolute value
    projEndLong = abs(dispVec * longVecE) / longVecMagE;
    projEndShort = abs(dispVec * shortVecE) / shortVecMagE;
    
    %get second short vector and project along it if problem is 3D
    if probDim == 3
        shortVecS3D = shortVecS3DAll(:,timeGap,iStart); %second short vectors
        shortVecE3D = shortVecE3DAll(:,timeGap,iEnd);
        shortVecMagS3D = sqrt(shortVecS3D' * shortVecS3D); %their magnitudes
        shortVecMagE3D = sqrt(shortVecE3D' * shortVecE3D);
        projStartShort3D = abs(dispVec * shortVecS3D) / shortVecMagS3D; %projection
        projEndShort3D = abs(dispVec * shortVecE3D) / shortVecMagE3D;
    else %if problem is 2D, make values zero
        shortVecMagS3D = 0;
        shortVecMagE3D = 0;
        projStartShort3D = 0;
        projEndShort3D = 0;
    end

    %calculate the vector connecting the centers of the two tracks
    cen2cenVec = trackCenter(iStart,:) - trackCenter(iEnd,:);
    cen2cenVecMag = sqrt(cen2cenVec * cen2cenVec');

    %decide whether this is a possible link based on the types of
    %the two tracks
    switch trackTypeE
        case 1 %if end is directed
            switch trackTypeS
                case 1 %if start is directed

                    %calculate the square sine of the angle between velocity vectors
                    sin2Angle = 1 - (longVecE' * longVecS / ...
                        (longVecMagE * longVecMagS))^2;
                    
                    %calculate the square sine of the angle between each
                    %motion direction vector and the center-to-center vector
                    sin2AngleE = 1 - (cen2cenVec * longVecE / ...
                        (longVecMagE * cen2cenVecMag))^2;
                    sin2AngleS = 1 - (cen2cenVec * longVecS / ...
                        (longVecMagS * cen2cenVecMag))^2;

                    %check whether the end of track iEnd is within the search
                    %rectangle of the start of track iStart and vice versa,
                    %whether the angle between the two directions of motion
                    %is within acceptable bounds, and whether the angle
                    %between directions of motion and vector connecting end
                    %and start is within acceptable bounds
                    possibleLink = ((projEndLong <= longVecMagE && ...
                        projEndShort <= shortVecMagE && ...
                        projEndShort3D <= shortVecMagE3D) || ...
                        (projStartLong <= longVecMagS && ...
                        projStartShort <= shortVecMagS && ...
                        projStartShort3D <= shortVecMagS3D)) && ...
                        sin2Angle <= sin2AngleMax && ...
                        (sin2AngleE <= sin2AngleMaxVD || ...
                        sin2AngleS <= sin2AngleMaxVD);
                    
                case 0 %if start is Brownian

                    %calculate the square sine of the angle between the
                    %end's motion direction vector and the displacement vector
                    sin2AngleE = 1 - (cen2cenVec * longVecE / ...
                        (longVecMagE * cen2cenVecMag))^2;

                    %check whether the start of track iStart is within the
                    %search rectangle of the end of track iEnd, whether the
                    %end of track iEnd is within the search disc of the start
                    %of track iStart, and whether the angle between directions 
                    %of motion and vector connecting end and start is 
                    %within acceptable bounds
                    possibleLink = ((projEndLong <= longVecMagE && ...
                        projEndShort <= shortVecMagE && ...
                        projEndShort3D <= shortVecMagE3D) || ...
                        dispVecMag <= longVecMagS) && ...
                        sin2AngleE <= sin2AngleMaxVD;

                otherwise %if start is undetermined

                    %calculate the square sine of the angle between the
                    %end's motion direction vector and the displacement vector
                    sin2AngleE = 1 - (cen2cenVec * longVecE / ...
                        (longVecMagE * cen2cenVecMag))^2;

                    %check whether the start of track iStart is within the search
                    %rectangle of the end of track iEnd, and whether the 
                    %angle between directions of motion and vector
                    %connecting end and start is within acceptable bounds
                    possibleLink = (projEndLong <= longVecMagE && ...
                        projEndShort <= shortVecMagE && ...
                        projEndShort3D <= shortVecMagE3D) && ...
                        sin2AngleE <= sin2AngleMaxVD;

            end
        case 0 %if end is Brownian
            switch trackTypeS
                case 1 %if start is directed

                    %calculate the square sine of the angle between the
                    %start's motion direction vector and the displacement vector
                    sin2AngleS = 1 - (cen2cenVec * longVecS / ...
                        (longVecMagS * cen2cenVecMag))^2;

                    %check whether the end of track iEnd is within the search
                    %rectangle of the start of track iStart, whether the
                    %start of track iStart is within the search disc of the
                    %end of track iEnd, and whether the angle between
                    %directions of motion and vector connecting end and
                    %start is within acceptable bounds
                    possibleLink = ((projStartLong <= longVecMagS && ...
                        projStartShort <= shortVecMagS && ...
                        projStartShort3D <= shortVecMagS3D) || ...
                        dispVecMag <= longVecMagE) && ...
                        sin2AngleS <= sin2AngleMaxVD;

                case 0 %if start is Brownian

                    %check whether the end of track iEnd is within the search
                    %disc of the start of track iStart and vice versa
                    possibleLink = (dispVecMag <= longVecMagE) || ...
                        (dispVecMag <= longVecMagS);

                otherwise %if start is undetermined

                    %check whether the start of track iStart is within the search
                    %disc of the end of track iEnd
                    possibleLink = dispVecMag <= longVecMagE;

            end
        otherwise %if end is undetermined

            switch trackTypeS
                case 1 %if start is directed

                    %calculate the square sine of the angle between the
                    %start's motion direction vector and the displacement vector
                    sin2AngleS = 1 - (cen2cenVec * longVecS / ...
                        (longVecMagS * cen2cenVecMag))^2;

                    %check whether the end of track iEnd is within the search
                    %rectangle of the start of track iStart
                    possibleLink = (projStartLong <= longVecMagS && ...
                        projStartShort <= shortVecMagS && ...
                        projStartShort3D <= shortVecMagS3D) && ...
                        sin2AngleS <= sin2AngleMaxVD;

                otherwise %if start is Brownian or undetermined

                    %check whether the end of track iEnd is within the search
                    %disc of the start of track iStart
                    possibleLink = dispVecMag <= longVecMagS;

            end
    end

    %if this is a possible link ...
    if possibleLink

        %calculate the cost of linking
        dispVecMag2 = dispVecMag ^ 2;
        if trackTypeE == 1 && trackTypeS == 1
            cost12 = dispVecMag2 * (1 + mean([sin2Angle sin2AngleE sin2AngleS])) ...
                / timeScalingLin(timeGap);
        elseif trackTypeE == 1
            cost12 = dispVecMag2 * (1 + sin2AngleE) / timeScalingLin(timeGap);
        elseif trackTypeS == 1
            cost12 = dispVecMag2 * (1 + sin2AngleS) / timeScalingLin(timeGap);
        else
            cost12 = dispVecMag2 / timeScalingBrown(timeGap);
        end

        %penalize cost for lifetime considerations
        if ~isempty(lftCdf)
            cost12 = cost12 / oneMinusLftCdf(trackEndTime(iStart)-trackStartTime(iEnd)+2);
        end

        %if the lifetime consideration does not make this link impossible
        if ~isinf(cost12)
            
            %penalize cost for gap length considerations
            cost12 = cost12 * gapPenalty^(timeGap-1);
            
            %add this cost to the list of costs
            cost(iPair) = cost12;

            %specify the location of this pair in the cost matrix
            indx1(iPair) = iEnd; %row number
            indx2(iPair) = iStart; %column number
            
        end

    end %(if possibleLink)

end %(for iPair = 1 : numPairs)

%keep only pairs that turned out to be possible
possiblePairs = find(indx1 ~= 0);
indx1 = indx1(possiblePairs);
indx2 = indx2(possiblePairs);
cost  = cost(possiblePairs);
clear possiblePairs

%% Merging and splitting

%define some merging and splitting variables
numMerge  =  0; %index counting merging events
indxMerge = []; %vector storing merging track number
altCostMerge = []; %vector storing alternative costs of not merging
numSplit  =  0; %index counting splitting events
indxSplit = []; %vector storing splitting track number
altCostSplit = []; %vector storing alternative costs of not splitting

%if merging and splitting are to be considered ...
if mergeSplit

    %get the absolute upper limit of acceptable displacements in one frame
    %as the maximum of (maximum velocity multiplied by probDim*linStdMult(1),
    %maximum maxSearchRadiusCG)
    maxDispAllowed = max(max(xyzVel(:)) * probDim * linStdMult(1), ...
        max(maxSearchRadius));

    %costs of merging

    %go over all track end times
    for endTime = 1 : numFrames-1

        %find tracks that end in this frame
        endsToConsider = tracksPerFrame(endTime).ends;
        
        %find tracks that start before or in this frame and end after this
        %frame
        mergesToConsider = intersect(vertcat(tracksPerFrame(1:endTime).starts),...
            vertcat(tracksPerFrame(endTime+1:end).ends));

        %get index indicating frame of merging
        timeIndx  = endTime*8;

        %calculate displacement between track ends and other tracks in the
        %next frame
        dispMat2 = createDistanceMatrix(coordEnd(endsToConsider,:), ...
            trackedFeatInfo(mergesToConsider,timeIndx+1:timeIndx+probDim));

        %find possible pairs
        [indxEnd2,indxMerge2] = find(dispMat2 <= maxDispAllowed);
        numPairs = length(indxEnd2);
        
        %clear memory
        clear dispMat2
        
        %map from indices to track indices
        indxEnd2 = endsToConsider(indxEnd2);
        indxMerge2 = mergesToConsider(indxMerge2);

        %reserve memory for cost vectors and related vectors
        indx1MS   = zeros(numPairs,1);
        indx2MS   = zeros(numPairs,1);
        costMS    = zeros(numPairs,1);
        altCostMS = zeros(numPairs,1);
        indxMSMS  = zeros(numPairs,1);
        
        %go over all possible pairs
        for iPair = 1 : numPairs

            %get indices of ending track and track it might merge with
            iEnd = indxEnd2(iPair);
            iMerge = indxMerge2(iPair);

            %determine the search ellipse of track iEnd
            longVecE = longVecEAll(:,1,iEnd);
            shortVecE = shortVecEAll(:,1,iEnd);

            %calculate the magnitudes of the long and short search vectors
            %of the end
            longVecMagE = sqrt(longVecE' * longVecE);
            shortVecMagE = sqrt(shortVecE' * shortVecE);

            %calculate the vector connecting the end of track iEnd to the
            %point of merging and compute its magnitude
            dispVec = coordEnd(iEnd,:) - trackedFeatInfo(iMerge,...
                timeIndx+1:timeIndx+probDim);
            dispVecMag = sqrt(dispVec * dispVec');
            
            %project the connecting vector onto the long and short vectors
            %of track iEnd and take absolute value
            projEndLong = abs(dispVec * longVecE) / longVecMagE;
            projEndShort = abs(dispVec * shortVecE) / shortVecMagE;

            %get the amplitude at the end of track iEnd
            ampE = ampEnd(iEnd);

            %get the amplitude of the merging track at the point of merging
            %and the point before it
            ampM = trackedFeatInfo(iMerge,8*endTime+4); %at point of merging
            ampM1 = trackedFeatInfo(iMerge,8*(endTime-1)+4); %just before merging
            
            %calculate the ratio of the amplitude after merging to the sum
            %of the amplitudes before merging
            ampRatio = ampM / (ampE + ampM1);
            
            %decide whether this is a possible link based on displacement,
            %directionality and amplitude ratio
            if trackType(iEnd) == 1 %if ending track is linear

                %get second short vector and project along it if problem is 3D
                if probDim == 3
                    shortVecE3D = shortVecE3DAll(:,1,iEnd);
                    shortVecMagE3D = sqrt(shortVecE3D' * shortVecE3D);
                    projEndShort3D = abs(dispVec * shortVecE3D) / shortVecMagE3D;
                else %if problem is 2D, make values zero
                    shortVecMagE3D = 0;
                    projEndShort3D = 0;
                end

                %calculate the vector connecting the centers of the two tracks
                cen2cenVec = trackCenter(iEnd,:) - trackCenter(iMerge,:);
                cen2cenVecMag = sqrt(cen2cenVec * cen2cenVec');

                %calculate the square sine of the angle between the
                %end's motion direction vector and the center-to-center vector
                sin2AngleE = 1 - (cen2cenVec * longVecE / ...
                    (longVecMagE * cen2cenVecMag))^2;

                %check whether the merging feature is within the search
                %region of the end of track iEnd, that the center-to-center
                %vector is reasonably well aligned with the directionality
                %of track iEnd, and that the amplitude ratio is
                %within acceptable limits
                possibleLink = projEndLong <= longVecMagE && ...
                    projEndShort <= shortVecMagE && ...
                    projEndShort3D <= shortVecMagE3D && ...
                    ampRatio >= minAmpRatio && ampRatio <= maxAmpRatio && ...
                    sin2AngleE <= sin2AngleMaxVD;

            else %if ending track is Brownian or undetermined

                %assign the dummy value of zero to sin2AngleE
                sin2AngleE = 0;
                
                %look at displacement and amplitude ratio only (no
                %directionality)
                possibleLink = dispVecMag <= longVecMagE && ...
                    ampRatio >= minAmpRatio && ampRatio <= maxAmpRatio;

            end

            %if this is a possible link ...
            if possibleLink

                %calculate the cost of linking
                dispVecMag2 = dispVecMag ^ 2; %due to displacement
                ampCost = ampRatio; %due to amplitude
                ampCost(ampCost<1) = ampCost(ampCost<1) ^ (-2); %punishment harsher when intensity of merged feature < sum of intensities of merging features
                cost12 = dispVecMag2 * ampCost * (1 + sin2AngleE); %cost

                %penalize cost for lifetime considerations
                if ~isempty(lftCdf)
                    cost12 = cost12 / oneMinusLftCdf(trackEndTime(iMerge)-trackStartTime(iEnd)+2);
                end

                %if the lifetime consideration does not make this link impossible
                if ~isinf(cost12)

                    %add this cost to the list of costs
                    costMS(iPair) = cost12;

                    %check whether the track being merged with has had
                    %something possibly merging with it in this same frame
                    prevAppearance = find(indxMSMS == iMerge);

                    %if this track in this frame did not appear before ...
                    if isempty(prevAppearance)

                        %increase the "merge index" by one
                        numMerge = numMerge + 1;

                        %save the merging track's index
                        indxMSMS(iPair) = iMerge;

                        %store the location of this pair in the cost matrix
                        indx1MS(iPair) = iEnd; %row number
                        indx2MS(iPair) = numMerge+numTracks; %column number

                        %calculate the alternative cost of not merging for the
                        %track that the end is possibly merging with

                        %get the average square displacement in this track
                        trackCoord = trackedFeatInfo(indxMSMS(iPair),:);
                        trackCoord = reshape(trackCoord',8,[]);
                        dispVecMag2 = (diff(trackCoord,1,2)).^2;
                        dispVecMag2 = nanmean(dispVecMag2,2);
                        dispVecMag2 = sum(dispVecMag2(1:2));

                        %calculate intensity cost if no merge happens
                        ampCost = ampM / ampM1;
                        ampCost(ampCost<1) = ampCost(ampCost<1) ^ (-2);

                        %calculate alternative cost
                        cost12 = dispVecMag2 * ampCost;

                        %add this cost to the list of alternative costs
                        altCostMS(iPair) = cost12;

                    else %if this track in this frame appeared before

                        %do not increase the "merge index" or save the merging
                        %track's index (they are already saved)

                        %store the location of this pair in the cost matrix
                        indx1MS(iPair) = iEnd; %row number
                        indx2MS(iPair) = indx2MS(prevAppearance); %column number

                        %no need to calculate and save the alternative cost
                        %since that is already saved from previous encounter

                    end %(if isempty(prevAppearance))
                    
                end %(if ~isinf(cost12))
                
            end %(if possibleLink)

        end %(for iPair = 1 : numPairs)

        %keep only pairs that turned out to be possible
        possiblePairs = find(indx1MS ~= 0);
        indx1MS   = indx1MS(possiblePairs);
        indx2MS   = indx2MS(possiblePairs);
        costMS    = costMS(possiblePairs);
        possibleMerges = find(indxMSMS ~= 0);
        indxMSMS  = indxMSMS(possibleMerges);
        altCostMS = altCostMS(possibleMerges);
        clear possiblePairs possibleMerges
        
        %append these vectors to overall cost vector and related vectors
        indx1 = [indx1; indx1MS];
        indx2 = [indx2; indx2MS];
        cost  = [cost; costMS];
        altCostMerge = [altCostMerge; altCostMS];
        indxMerge = [indxMerge; indxMSMS];

    end %(for endTime = 1 : numFrames-1)

    %costs of splitting

    %go over all track starting times
    for startTime = 2 : numFrames

        %find tracks that start in this frame
        startsToConsider = tracksPerFrame(startTime).starts;
        
        %find tracks that start before this frame and end after or in this frame
        splitsToConsider = intersect(vertcat(tracksPerFrame(1:startTime-1).starts),...
            vertcat(tracksPerFrame(startTime:end).ends));

        %get index indicating time of splitting
        timeIndx  = (startTime-2)*8;

        %calculate displacement between track starts and other tracks in the
        %previous frame
        dispMat2 = createDistanceMatrix(coordStart(startsToConsider,:), ...
            trackedFeatInfo(splitsToConsider,timeIndx+1:timeIndx+probDim));

        %find possible pairs
        [indxStart2,indxSplit2] = find(dispMat2 <= maxDispAllowed);
        numPairs = length(indxStart2);
        
        %clear memory
        clear dispMat2

        %map from indices to track indices
        indxStart2 = startsToConsider(indxStart2);
        indxSplit2 = splitsToConsider(indxSplit2);

        %reserve memory for cost vectors and related vectors
        indx1MS   = zeros(numPairs,1);
        indx2MS   = zeros(numPairs,1);
        costMS    = zeros(numPairs,1);
        altCostMS = zeros(numPairs,1);
        indxMSMS  = zeros(numPairs,1);
        
        %go over all possible pairs
        for iPair = 1 : numPairs

            %get indices of starting track and track it might have split from
            iStart = indxStart2(iPair);
            iSplit = indxSplit2(iPair);

            %determine the search ellipse of track iStart
            longVecS = longVecSAll(:,1,iStart);
            shortVecS = shortVecSAll(:,1,iStart);

            %calculate the magnitudes of the long and short search vectors
            %of the start
            longVecMagS = sqrt(longVecS' * longVecS);
            shortVecMagS = sqrt(shortVecS' * shortVecS);

            %calculate the vector connecting the end of track iStart to the
            %point of splitting and compute its magnitude
            dispVec = coordStart(iStart,:) - trackedFeatInfo(iSplit,...
                timeIndx+1:timeIndx+probDim);
            dispVecMag = sqrt(dispVec * dispVec');

            %project the connecting vector onto the long and short vectors
            %of track iStart and take absolute value
            projStartLong = abs(dispVec * longVecS) / longVecMagS;
            projStartShort = abs(dispVec * shortVecS) / shortVecMagS;

            %get the amplitude at the start of track iStart
            ampS = ampStart(iStart);

            %get the amplitude of the splitting track at the point of splitting
            %and the point before it
            ampSp1 = trackedFeatInfo(iSplit,8*(startTime-1)+4); %at point of splitting
            ampSp = trackedFeatInfo(iSplit,8*(startTime-2)+4); %just before splitting
            
            %calculate the ratio of the amplitude before splitting to the sum
            %of the amplitudes after splitting
            ampRatio = ampSp / (ampS + ampSp1);
            
            %decide whether this is a possible link based on displacement,
            %directionality and amplitude ratio
            if trackType(iStart) == 1

                %get second short vector and project along it if problem is 3D
                if probDim == 3
                    shortVecS3D = shortVecS3DAll(:,1,iStart); %second short vectors
                    shortVecMagS3D = sqrt(shortVecS3D' * shortVecS3D); %their magnitudes
                    projStartShort3D = abs(dispVec * shortVecS3D) / shortVecMagS3D; %projection
                else %if problem is 2D, make values zero
                    shortVecMagS3D = 0;
                    projStartShort3D = 0;
                end

                %calculate the vector connecting the centers of the two tracks
                cen2cenVec = trackCenter(iStart,:) - trackCenter(iSplit,:);
                cen2cenVecMag = sqrt(cen2cenVec * cen2cenVec');

                %calculate the square sine of the angle between the
                %end's motion direction vector and the center-to-center vector
                sin2AngleS = 1 - (cen2cenVec * longVecS / ...
                    (longVecMagS * cen2cenVecMag))^2;

                %check whether the merging feature is within the search
                %region of the end of track iEnd, that the center-to-center
                %vector is reasonably well aligned with the directionality
                %of track iEnd, and that the amplitude ratio is
                %within acceptable limits
                possibleLink = projStartLong <= longVecMagS && ...
                    projStartShort <= shortVecMagS && ...
                    projStartShort3D <= shortVecMagS3D && ...
                    ampRatio >= minAmpRatio && ampRatio <= maxAmpRatio && ...
                    sin2AngleS <= sin2AngleMaxVD;
            else
                
                %assign the dummy value of zero to sin2AngleS
                sin2AngleS = 0;
                
                %look at displacement and amplitude ratio only (no
                %directionality)
                possibleLink = dispVecMag <= longVecMagS && ...
                    ampRatio >= minAmpRatio && ampRatio <= maxAmpRatio;
                
            end

            %if this is a possible link ...
            if possibleLink

                %calculate the cost of linking
                dispVecMag2 = dispVecMag ^ 2; %due to displacement
                ampCost = ampRatio; %due to amplitude
                ampCost(ampCost<1) = ampCost(ampCost<1) ^ (-2); %punishment harsher when intensity of splitting feature < sum of intensities of features after splitting
                cost12 = dispVecMag2 * ampCost * (1 + sin2AngleS);

                %penalize cost for lifetime considerations
                if ~isempty(lftCdf)
                    cost12 = cost12 / oneMinusLftCdf(trackEndTime(iStart)-trackStartTime(iSplit)+2);
                end

                %if the lifetime consideration does not make this link impossible
                if ~isinf(cost12)

                    %add this cost to the list of costs
                    costMS(iPair) = cost12;

                    %check whether the track being split from has had something
                    %possibly splitting from it in this same frame
                    prevAppearance = find(indxMSMS == iSplit);

                    %if this track in this frame did not appear before ...
                    if isempty(prevAppearance)

                        %increase the "split index" by one
                        numSplit = numSplit + 1;

                        %save the splitting track's number
                        indxMSMS(iPair) = iSplit;

                        %store the location of this pair in the cost matrix
                        indx1MS(iPair) = numSplit+numTracks; %row number
                        indx2MS(iPair) = iStart; %column number

                        %calculate the alternative cost of not splitting for the
                        %track that the start is possibly splitting from

                        %get the average square displacement in this track
                        trackCoord = trackedFeatInfo(indxMSMS(iPair),:);
                        trackCoord = reshape(trackCoord',8,[]);
                        dispVecMag2 = (diff(trackCoord,1,2)).^2;
                        dispVecMag2 = nanmean(dispVecMag2,2);
                        dispVecMag2 = sum(dispVecMag2(1:2));

                        %calculate intensity cost if no split happens
                        ampCost = ampSp / ampSp1;
                        ampCost(ampCost<1) = ampCost(ampCost<1) ^ (-2);

                        %calculate alternative cost
                        cost12 = dispVecMag2 * ampCost;

                        %add this cost to the list of alternative costs
                        altCostMS(iPair) = cost12;

                    else %if this track in this frame appeared before

                        %do not increase the "split index" or save the
                        %splitting track's index (they are already saved)

                        %store the location of this pair in the cost matrix
                        indx1MS(iPair) = indx1MS(prevAppearance); %row number
                        indx2MS(iPair) = iStart; %column number

                        %no need to calculate and save the alternative cost
                        %since that is already saved from previous appearance

                    end %(if isempty(prevAppearance))
                    
                end %(if ~isinf(cost12))

            end %(if possibleLink)

        end %(for for iPair = 1 : numPairs)

        %keep only pairs that turned out to be possible
        possiblePairs = find(indx1MS ~= 0);
        indx1MS   = indx1MS(possiblePairs);
        indx2MS   = indx2MS(possiblePairs);
        costMS    = costMS(possiblePairs);
        possibleSplits = find(indxMSMS ~= 0);
        altCostMS = altCostMS(possibleSplits);
        indxMSMS  = indxMSMS(possibleSplits);
        clear possiblePairs possibleSplits
        
        %append these vectors to overall cost and related vectors
        indx1 = [indx1; indx1MS];
        indx2 = [indx2; indx2MS];
        cost  = [cost; costMS];
        altCostSplit = [altCostSplit; altCostMS];
        indxSplit = [indxSplit; indxMSMS];

    end %(for startTime = 2 : numFrames)

end %(if mergeSplit)

%create cost matrix without births and deaths
numEndSplit = numTracks + numSplit;
numStartMerge = numTracks + numMerge;
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
