function tracksCrop = cropCompTracksFR(tracks0,frameRange)
%cropCompTracksFR crops tracks to a particular frame range
%
%SYNOPSIS tracksCrop = cropCompTracksFR(tracks0,frameRange)
%
%INPUT  tracks0      : Original tracks, in form of output of
%                      trackCloseGapsKalman.
%       frameRange   : Row vector with two entries indicating time range to
%                      retain.
%
%OUTPUT tracksCrop   : Same as input tracks, just cropped in time.
%
%REMARKS Code still needs to handle case when a compound track gets broken
%        into multiple non-interacting tracks.
%
%Khuloud Jaqaman, December 2014

%% Input

if nargin < 2
    error('cropCompTracksFR: Incorrect number of input arguments!')
end

%determine which frames to keep
tpKeep = frameRange(1) : frameRange(2);

%thus determine which columns to keep in tracksCoordAmpCG
colKeep = (repmat(tpKeep,8,1)-1)*8 + repmat((1:8)',1,length(tpKeep));
colKeep = colKeep(:);

%check if tracks are in sparse form
sparseForm = issparse(tracks0(1).tracksFeatIndxCG);

%% Cropping

%get number of tracks
numTracks = length(tracks0);

%go over all tracks and crop frames
tracksCrop = tracks0;
for iTrack = 1 : numTracks
    
    %convert the current track's information to matrix format
    [trackedFeatureInfo,trackedFeatureIndx] = convStruct2MatIgnoreMS(tracks0(iTrack));
    
    %keep only the time points of interest
    trackedFeatureIndx = trackedFeatureIndx(:,tpKeep);
    trackedFeatureInfo = trackedFeatureInfo(:,colKeep);

    %convert zeros to NaNs if original tracks were sparse
    if sparseForm
        trackedFeatureInfo(trackedFeatureInfo==0) = NaN;
    end
    
    %get each track segment's new start, end and life time
    segSEL = getTrackSEL(trackedFeatureInfo);
    numSeg = size(segSEL,1);
    
    %find segments that survive the cropping and those that do not
    indxStay = find(~isnan(segSEL(:,3)));
    indxGone = setdiff((1:numSeg)',indxStay);
    if isempty(indxGone)
        indxGone  = [];
    end
    
    %get the track's original sequence of events
    seqOfEvents = tracks0(iTrack).seqOfEvents;
    
    %"delete" merges and splits happening outside of the frame range
    tmp = ~isnan(seqOfEvents(:,4)) & ( seqOfEvents(:,1)<frameRange(1) | seqOfEvents(:,1)>frameRange(2) );
    seqOfEvents(tmp,4) = NaN;
    
    %assign new start and end times - merges will need an addition of one,
    %done two steps down
    for iSeg = 1 : numSeg
        rowsSeg = find(seqOfEvents(:,3)==iSeg);
        seqOfEvents(rowsSeg(1),1) = segSEL(iSeg,1);
        seqOfEvents(rowsSeg(2),1) = segSEL(iSeg,2);
    end
    
    %replace segments that did not survive the cropping with NaN in both
    %the 3rd and 4th column
    for iSeg = indxGone'
        seqOfEvents(seqOfEvents(:,3)==iSeg,3) = NaN;
        seqOfEvents(seqOfEvents(:,4)==iSeg,4) = NaN;
    end
    
    %for surviving merges, add one to end time to follow convention
    rowsSeg = find(seqOfEvents(:,2)==2 & ~isnan(seqOfEvents(:,4)));
    seqOfEvents(rowsSeg,1) = seqOfEvents(rowsSeg,1) + 1;
    
    %keep only rows that belong to surviving segments
    seqOfEvents = seqOfEvents(~isnan(seqOfEvents(:,3)),:);
    trackedFeatureInfo = trackedFeatureInfo(indxStay,:);
    trackedFeatureIndx = trackedFeatureIndx(indxStay,:);
    
    %renumber segments to reflect new trackedFeatureInfo and
    %trackedFeatureIndx
    for iStay = 1 : length(indxStay)
        iSeg = indxStay(iStay);
        seqOfEvents(seqOfEvents(:,3)==iSeg,3) = iStay;
        seqOfEvents(seqOfEvents(:,4)==iSeg,4) = iStay;
    end
    
    %convert to sparse if input was sparse
    if sparseForm
        trackedFeatureIndx = sparse(trackedFeatureIndx);
        trackedFeatureInfo(isnan(trackedFeatureInfo)) = 0;
        trackedFeatureInfo = sparse(trackedFeatureInfo);
    end
    
    %store the cropped compound track
    tracksCrop(iTrack).tracksFeatIndxCG = trackedFeatureIndx;
    tracksCrop(iTrack).tracksCoordAmpCG = trackedFeatureInfo;
    tracksCrop(iTrack).seqOfEvents = seqOfEvents;
    
end

%% ~~~ the end ~~~

