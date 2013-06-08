function numNeighbors = numNeighborsTrack(tracksFinal,winFrames,nRadius)
%NUMNEIGHBORSTRACK calculates number of neighbors for a molecule within a certain area
%
%SYNOPSIS numNeighbors = numNeighborsTrack(tracksFinal,nRadius)
%
%INPUT  tracksFinal : Output of trackCloseGapsKalman.
%       winFrames   : The SPT frames at which there are windows.
%       nRadius     : Radius to count number of neighbors.
%
%OUTPUT numNeighbors: Structure array with the field "value" storing number
%                     of neighbors per track segment.
%
%Khuloud Jaqaman, June 2013

%% Input

if nargin < 2
    disp('--numNeighborsTrack: Incorrect number of input arguments!');
    return
end

if nargin < 3 || isempty(nRadius)
    nRadius = 5;
end

%generate winFrameMin and winFramesMax for for loop
winFramesMid = floor((winFrames(1:end-1) + winFrames(2:end))/2);
winFramesMin = [1 winFramesMid];
winFramesMax = [winFramesMid winFrames(end)+1];
numWinFrames = length(winFramesMin);

%% Number of neighbors calculation

%get number of compound tracks
numTracks = length(tracksFinal);

%initilize structure arrays
[meanCoord,meanTime,numNeighbors] = deal(repmat(struct('value',[]),numTracks,1));

%get each track segment's mean position and time 
for iTrack = 1 : numTracks
    xCoordAll = tracksFinal(iTrack).tracksCoordAmpCG(:,1:8:end);
    yCoordAll = tracksFinal(iTrack).tracksCoordAmpCG(:,2:8:end);
    trackSEL = getTrackSEL(tracksFinal(iTrack),1);
    meanCoord(iTrack).value = [nanmean(xCoordAll,2) nanmean(yCoordAll,2)];
    meanTime(iTrack).value = nanmean(trackSEL(:,1:2),2);
end
meanCoordAll = vertcat(meanCoord.value);
meanTimeAll = vertcat(meanTime.value);

%group the track segments based on what window frame range they fall into
trackGroup = repmat(struct('indx',[]),numWinFrames,1);
for iWinFrame = 1 : numWinFrames
  
    %get current spt frame number and next spt frame number
    minFrame = winFramesMin(iWinFrame);
    maxFrame = winFramesMax(iWinFrame);
        
    %find track segments whose "average" time is between minFrame and maxFrame
    indxFrameRange = find(meanTimeAll>=minFrame & meanTimeAll<maxFrame);

    %store this information for later use
    trackGroup(iWinFrame).indx = indxFrameRange;
    
end

%calculate number of neighbors within nRadius
for iTrack = 1 : numTracks
    for iSeg = 1 : size(meanCoord(iTrack).value);
        
        %get current track segment's information
        trackCoord = meanCoord(iTrack).value(iSeg,:);
        trackTime = meanTime(iTrack).value(iSeg);
        
        %determine what tmie group it belongs to and get concurrent track
        %segments
        groupIndx = find(trackTime<winFramesMax,1,'first');
        trackIndx = trackGroup(groupIndx).indx;
        
        %calculate distance between track segments
        distMat = createDistanceMatrix(trackCoord,meanCoordAll(trackIndx,:));
        
        %get number of neighbors within nRadius
        tmpNN = length(find(distMat < nRadius)) - 1;
        if groupIndx==1 || groupIndx == numWinFrames
            tmpNN = tmpNN * 2;
        end
        numNeighbors(iTrack).value(iSeg,1) = tmpNN;
        
    end
end

%% ~~~ the end ~~~

