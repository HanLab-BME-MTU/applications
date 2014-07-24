function numNeighbors = numNeighborsTrack(tracksFinal,spatialRadius,timeRadius)
%NUMNEIGHBORSTRACK calculates number of neighbors for a molecule within a certain area
%
%SYNOPSIS numNeighbors = numNeighborsTrack(tracksFinal,spatialRadius)
%
%INPUT  tracksFinal  : Output of trackCloseGapsKalman.
%       spatialRadius: Spatial sadius to count number of neighbors.
%                      Optional. Default: 5 [whatever units].
%       timeRadius   : Time radius to count number of neighbors.
%                      Optional. Default: 5 frames.
%
%OUTPUT numNeighbors: Structure array with the field "value" storing number
%                     of neighbors per track segment.
%
%Khuloud Jaqaman, June 2013

%% Input

if nargin < 1
    disp('--numNeighborsTrack: Incorrect number of input arguments!');
    return
end

if nargin < 2 || isempty(spatialRadius)
    spatialRadius = 5;
end

if nargin < 3 || isempty(timeRadius)
    timeRadius = 5;
end

%get number of frames in movie
seqOfEvents = vertcat(tracksFinal.seqOfEvents);
numFrames = max(seqOfEvents(:,1));

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

%calculate number of neighbors within spatialRadius and timeRadius
for iTrack = 1 : numTracks
    for iSeg = 1 : size(meanCoord(iTrack).value);
        
        %get current track segment's information
        trackCoord = meanCoord(iTrack).value(iSeg,:);
        trackTime = meanTime(iTrack).value(iSeg);
        
        %determine the track segments that fall within timeRadius
        trackIndx = find( abs(meanTimeAll-trackTime) <= timeRadius );
        
        %calculate the distance between track segments
        distMat = createDistanceMatrix(trackCoord,meanCoordAll(trackIndx,:)); %#ok<FNDSB>
        
        %get number of neighbors within spatialRadius
        tmpNN = length(find(distMat <= spatialRadius)) - 1;
        
        %compensate in case track segment is too close to movie start/end
        timeFromStartEnd = min(trackTime-1,numFrames-trackTime);
        if timeFromStartEnd < timeRadius
            tmpNN = tmpNN * (2-timeFromStartEnd/timeRadius);
        end
        numNeighbors(iTrack).value(iSeg,1) = tmpNN;
        
    end
end

%% ~~~ the end ~~~

