function [costMat,nonlinkMarker,errFlag] = costMatCloseGapsSR(trackedFeatInfo,...
    trackStartTime,trackEndTime,searchRadius,timeWindow,gapPenalty,probDim)
%costMatCloseGapsSR provides a cost matrix for closing gaps in the tracks of static features for super-resolution applications
%
%SYNOPSIS [costMat,nonlinkMarker,errFlag] = costMatCloseGapsSR(trackedFeatInfo,...
%    trackedFeatIndx,trackStartTime,trackEndTime,searchRadius,timeWindow,...
%    gapPenalty,probDim)
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
%       trackStartTime : Starting time of all tracks.
%       trackEndTime   : Ending time of all tracks.
%       searchRadius   : Search radius for closing gaps, in pixels.
%       timeWindow     : Largest time gap between the end of a track and
%                        the beginning of another that could be connected
%                        to it.
%       gapPenalty     : Penalty for increasing temporary disappearance
%                        time, to be used in gap closing cost. Disappearing
%                        for n frames, i.e. closing a gap of n+1 frames,
%                        gets a penalty of gapPenalty^n.
%       probDim        : Problem dimensionality. 2 (for 2D) or 3 (for 3D).
%
%OUTPUT costMat       : Cost matrix.
%       nonlinkMarker : Value indicating that a link is not allowed.
%       errFlag       : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, August 2011

%% Output

costMat = [];
nonlinkMarker = [];
errFlag = [];

%% Input

%check whether correct number of input arguments was used
if nargin ~= nargin('costMatCloseGapsSR')
    disp('--costMatCloseGapsSR: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

%find the number of tracks to be linked and the number of frames in the movie
[numTracks,numFrames] = size(trackedFeatInfo);
numFrames = numFrames / 8;

%list the tracks that start and end in each frame
tracksPerFrame = repmat(struct('starts',[],'ends',[]),numFrames,1);
for iFrame = 1 : numFrames    
    tracksPerFrame(iFrame).starts = find(trackStartTime == iFrame); %starts
    tracksPerFrame(iFrame).ends = find(trackEndTime == iFrame); %ends
end

%% Pre-processing

%get the x,y-coordinates and amplitudes at the starts of tracks
coordStart = zeros(numTracks,probDim);
ampStart   = zeros(numTracks,1);
for iTrack = 1 : numTracks
    coordStart(iTrack,:) = full(trackedFeatInfo(iTrack,...
        (trackStartTime(iTrack)-1)*8+1:(trackStartTime(iTrack)-1)*8+probDim));
    ampStart(iTrack) = full(trackedFeatInfo(iTrack,(trackStartTime(iTrack)-1)*8+4));
end

%get the x,y-coordinates and amplitudes at the ends of tracks
coordEnd = zeros(numTracks,probDim);
ampEnd   = zeros(numTracks,1);
for iTrack = 1 : numTracks
    coordEnd(iTrack,:) = full(trackedFeatInfo(iTrack,...
        (trackEndTime(iTrack)-1)*8+1:(trackEndTime(iTrack)-1)*8+probDim));
    ampEnd(iTrack) = full(trackedFeatInfo(iTrack,(trackEndTime(iTrack)-1)*8+4));
end

%get the average coordinate of each track
coordMean = zeros(numTracks,probDim);
for iTrack = 1 : numTracks
    xCoordTmp = full(trackedFeatInfo(iTrack,...
        (trackStartTime(iTrack)-1)*8+1:8:(trackEndTime(iTrack)-1)*8+1));
    yCoordTmp = full(trackedFeatInfo(iTrack,...
        (trackStartTime(iTrack)-1)*8+2:8:(trackEndTime(iTrack)-1)*8+2));
    zCoordTmp = full(trackedFeatInfo(iTrack,...
        (trackStartTime(iTrack)-1)*8+3:8:(trackEndTime(iTrack)-1)*8+3));
    xyzCoordTmp = [mean(xCoordTmp,2) mean(yCoordTmp,2) mean(zCoordTmp,2)];
    coordMean(iTrack,:) = xyzCoordTmp(1:probDim);
end

%% Gap closing

%find all pairs of ends and starts that can potentially be linked
%determine this by looking at time gaps between ends and starts
%and by looking at the distance between pairs
indxEnd2 = [];
indxStart2 = [];

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
        %         dispMat2 = createDistanceMatrix(coordMean(endsToConsider,:),...
        %             coordMean(startsToConsider,:));
        
        %find possible pairs
        [indxEnd3,indxStart3] = find(dispMat2 <= searchRadius);
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
clear dispMat2

%reserve memory for cost matrix vectors
indx1 = zeros(numPairs,1); %row number in cost matrix
indx2 = zeros(numPairs,1); %column number in cost matrix
cost  = zeros(numPairs,1); %cost value

%go over all possible pairs of starts and ends
for iPair = 1 : numPairs
    
    %get indices of starts and ends
    iStart = indxStart2(iPair);
    iEnd = indxEnd2(iPair);
    
    %determine the time gap between them
    timeGap = trackStartTime(iStart) - trackEndTime(iEnd);
    
    %calculate the vector connecting the end of track iEnd to the
    %start of track iStart and compute its magnitude
    dispVec = coordEnd(iEnd,:) - coordStart(iStart,:);
    %     dispVec = coordMean(iEnd,:) - coordMean(iStart,:);
    dispVecMag = norm(dispVec);
    
    %check whether the end of track iEnd is within the search
    %disc of the start of track iStart
    possibleLink = dispVecMag <= searchRadius;
    
    %if this is a possible link ...
    if possibleLink
        
        cost12 = dispVecMag^2;
        
        %penalize cost for gap length considerations
        cost12 = cost12 * gapPenalty^(timeGap-1);
        
        %add this cost to the list of costs
        cost(iPair) = cost12;
        
        %specify the location of this pair in the cost matrix
        indx1(iPair) = iEnd; %row number
        indx2(iPair) = iStart; %column number
        
    end %(if possibleLink)
    
end %(for iPair = 1 : numPairs)

%keep only pairs that turned out to be possible
possiblePairs = find(indx1 ~= 0);
indx1 = indx1(possiblePairs);
indx2 = indx2(possiblePairs);
cost  = cost(possiblePairs);
costMat = sparse(indx1,indx2,cost,numTracks,numTracks);

clear possiblePairs

%% Append cost matrix to allow births and deaths ...

%determine the cost of birth and death
costBD = 2*max(max(costMat(:)),eps);

%get the cost for the lower right block
costLR = min(min(costMat(:))-1,-1);

costMat = [costMat ... %costs for links (gap closing + merge/split)
    spdiags(costBD*ones(numTracks,1),0,numTracks,numTracks); ... %costs for death
    spdiags(costBD*ones(numTracks,1),0,numTracks,numTracks) ...  %costs for birth
    sparse(indx2,indx1,costLR*ones(length(indx1),1),numTracks,numTracks)]; %dummy costs to complete the cost matrix

%determine the nonlinkMarker
nonlinkMarker = min(floor(full(min(min(costMat))))-5,-5);


%% ~~~ the end ~~~
