function [movieInfo,tracksSimMiss,percentMissingGT] = genMovieInfoFromMtTracks(tracksSim,percentMissing)
%GENMOVIEINFOFROMTRACKS2 generates a list of detected features per frame from supplied tracks
%
%SYNOPSIS [movieInfo,tracksSimMiss] = genMovieInfoFromMtTracks(tracksSim,percentMissing)
%
%INPUT  tracksSim     : Output of simulateMimickCD36_MS.
%       percentMissing: Percentage of missing features in movie.
%
%OUTPUT movieInfo: List of detected features per frame, in the format
%                  required for the input of trackWithGapClosing and
%                  trackCloseGapsKalman.
%
%Khuloud Jaqaman, October 2007
% Modified by Sebastien Beson

%get number of frames in movie
seqOfEvents = vertcat(tracksSim.seqOfEvents);
numFrames = max(seqOfEvents(:,1));

%get number of tracks
numTracks = length(tracksSim);

%define standard deviation of missing features
missStd = 0.1*percentMissing;

%pre-allocate memory for movieInfo
movieInfo = repmat(struct('xCoord',[],'yCoord',[],'amp',[]),numFrames,1);

%get number of segments making each track
numSegments = zeros(numTracks,1);
for iTrack = 1 : numTracks
    numSegments(iTrack) = size(tracksSim(iTrack).tracksCoordAmpCG,1);
end

%locate the row of the first segment of each compound track in the
%big matrix of all tracks (to be constructed in the next step)
trackStartRow = ones(numTracks,1);
for iTrack = 2 : numTracks
    trackStartRow(iTrack) = trackStartRow(iTrack-1) + numSegments(iTrack-1);
end

%% put all tracks together in a big matrix
%store coordinates just before/after a merge/split as their negative values
%in order to distinguish them from coordinates not related to merges/splits
%also store in negative values the coordinates at the start/end of a
%segment
trackedFeatureInfo = NaN*ones(trackStartRow(end)+numSegments(end)-1,8*numFrames);

for iTrack = 1 : numTracks

    %get track's sequence of events and coordinates
    seqOfEvents = tracksSim(iTrack).seqOfEvents;
    tracksCoordAmpCG = tracksSim(iTrack).tracksCoordAmpCG;

    %get track's start time and end time
    startTime = seqOfEvents(1,1);
    endTime = seqOfEvents(end,1);

    %find starts and ends of segments making this compound track
    indxSE = find(isnan(seqOfEvents(:,4)))';
    
    %go over all starts/ends
    for iSE = indxSE
        
        %get frame of start/end and track segment involved
        frameSE = seqOfEvents(iSE,1);
        segment1 = seqOfEvents(iSE,3);
        
        %make the coordinates at the start/end negative
        tracksCoordAmpCG(segment1,(frameSE-startTime)*8+1:(frameSE-startTime+1)*8) = ...
            -abs(tracksCoordAmpCG(segment1,(frameSE-startTime)*8+1:(frameSE-startTime+1)*8));
        
    end
        
    %find merges and splits happening in this compound track
    indxMS = find(~isnan(seqOfEvents(:,4)))';

    %go over all merges/splits
    for iMS = indxMS

        %get frame of merge/split and track segments involved
        frameMS = seqOfEvents(iMS,1);
        segment1 = seqOfEvents(iMS,3);
        segment2 = seqOfEvents(iMS,4);

        %make the coordinates just before/after a merge/split negative
        tracksCoordAmpCG(segment1,(frameMS-startTime-1)*8+1:(frameMS-startTime+1)*8) = ...
            -abs(tracksCoordAmpCG(segment1,(frameMS-startTime-1)*8+1:(frameMS-startTime+1)*8));
        tracksCoordAmpCG(segment2,(frameMS-startTime-1)*8+1:(frameMS-startTime+1)*8) = ...
            -abs(tracksCoordAmpCG(segment2,(frameMS-startTime-1)*8+1:(frameMS-startTime+1)*8));

    end

    %store the compound track's coordinates in big matrix
    trackedFeatureInfo(trackStartRow(iTrack):trackStartRow(iTrack)+...
        numSegments(iTrack)-1,8*(startTime-1)+1:8*endTime) = tracksCoordAmpCG;

end

%% in each frame, delete "percentMissing+-std" of the features that can be deleted and store information in movieInfo
%(i.e. features not just before/after a merge/split)

trackedFeatureInfoMiss = trackedFeatureInfo;
nGaps = 0;
%go over all frames ...
for iFrame = 1 : numFrames

    %find number of features in this frame
    numFeat = length(find(~isnan(trackedFeatureInfoMiss(:,(iFrame-1)*8+1))));
    
    %find indices of features that can be deleted in this frame
    indxCanDelete = find(trackedFeatureInfoMiss(:,(iFrame-1)*8+1)>0);
    
    %determine number of features to delete
    numFeatDelete = round((percentMissing+randn(1)*missStd)*numFeat/100);

    %randomly choose which features to delete
    if length(indxCanDelete) > numFeatDelete
        deleteIndx = randsample(indxCanDelete,numFeatDelete);
    else
        deleteIndx = indxCanDelete;
    end

    %delete features by making their coordinates NaN
    trackedFeatureInfoMiss(deleteIndx,(iFrame-1)*8+1:iFrame*8) = NaN;
    nGaps = nGaps + length(deleteIndx);
end

%% check for isolated features
allFrames = 1:numFrames;
trackedFeatureInfoMiss2 =trackedFeatureInfoMiss;
for iTrack =1:size(trackedFeatureInfoMiss)
    
    features =~isnan(trackedFeatureInfoMiss(iTrack,(0:numFrames-1)*8+1));
    positions = 2:numFrames-1;
    isIsolatedFeature = @(x) ~features(x-1) && features(x) && ~features(x+1);
    isolatedFeatPos =positions(arrayfun(isIsolatedFeature,positions));
    for i = isolatedFeatPos
        trackedFeatureInfoMiss(iTrack,(i-2)*8+1:(i-1)*8)= ...
            trackedFeatureInfo(iTrack,(i-2)*8+1:(i-1)*8);
         trackedFeatureInfoMiss(iTrack,(i)*8+1:(i+1)*8)= ...
            trackedFeatureInfo(iTrack,(i)*8+1:(i+1)*8);       
    end
end
for iTrack =1:size(trackedFeatureInfoMiss)
    
    features =~isnan(trackedFeatureInfoMiss(iTrack,(0:numFrames-1)*8+1));
    positions = 2:numFrames-1;
    isIsolatedFeature = @(x) ~features(x-1) && features(x) && ~features(x+1);
    isolatedFeatPos =positions(arrayfun(isIsolatedFeature,positions));
    assert(isempty(isolatedFeatPos));
end
percentMissingGT=nGaps/numFrames;
%% construct movieInfo

%go over all frames ...
for iFrame = 1 : numFrames
    %get the coordinates of features detected in this frame
    coord = trackedFeatureInfoMiss(:,(iFrame-1)*8+1:(iFrame-1)*8+4);
    coord = abs(coord(~isnan(coord(:,1)),:));
    numFeat = size(coord,1);

    %store feature information in movieInfo
    movieInfo(iFrame).xCoord = [coord(:,1) zeros(numFeat,1)];
    movieInfo(iFrame).yCoord = [coord(:,2) zeros(numFeat,1)];
    movieInfo(iFrame).amp    = [coord(:,4) zeros(numFeat,1)];

end


%% construct tracksSimMiss which includes NaNs for "closed gaps"

%initialize tracksSimMiss
tracksSimMiss = tracksSim;

for iTrack = 1 : numTracks
    
    %get track's sequence of events
    seqOfEvents = tracksSimMiss(iTrack).seqOfEvents;

    %get track's start time and end time
    startTime = seqOfEvents(1,1);
    endTime = seqOfEvents(end,1);

    %extract track's coordinates from big matrix
    tracksCoordAmpCG = trackedFeatureInfoMiss(trackStartRow(iTrack):...
        trackStartRow(iTrack)+numSegments(iTrack)-1,(startTime-1)*8+1:...
        endTime*8);
    tracksCoordAmpCG = abs(tracksCoordAmpCG);
    
    %store coordinates in tracksSimMiss
    tracksSimMiss(iTrack).tracksCoordAmpCG = tracksCoordAmpCG;
    
end


%% %%% ~~ the end ~~ %%%
