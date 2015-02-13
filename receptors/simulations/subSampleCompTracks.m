function tracksSub = subSampleCompTracks(tracks0,timeStep0,timeStepSub)
%SUBSAMPLECOMPTRACKS sub-samples compound tracks to a requested sampling rate
%
%SYNOPSIS tracksSub = subSampleCompTracks(tracks0,timeStep0,timeStepSub)
%
%INPUT  tracks0      : Original tracks, in form of output of
%                      trackCloseGapsKalman.
%       timeStep0    : Time step in original tracks.
%       timeStepSub  : Time step to sample down to, preferably an integer
%                      multiple of timeStep0. If not, it will get
%                      approximated to the closest integer.
%
%OUTPUT tracksSub    : Same as input tracks, just sub-sampled in time.
%
%REMARKS Code still needs to handle case when a compound track gets broken
%        into multiple non-interacting tracks.
%
%Khuloud Jaqaman, February 2013

%% Input

if nargin < 3
    error('subSampleCompTracks: Incorrect number of input arguments!')
end

if timeStepSub < timeStep0
    error('subSampleCompTracks: timeStepSub should be larger than timeStep0!')
end

%calculate sub-sampling factor
subSampFact = round(timeStepSub/timeStep0);

%determine which time points to keep
seqOfEvents = vertcat(tracks0.seqOfEvents);
numTP0 = max(seqOfEvents(:,1));
tpKeep = 1 : subSampFact : numTP0;

%thus determine which columns to keep in tracksCoordAmpCG
colKeep = (repmat(tpKeep,8,1)-1)*8 + repmat((1:8)',1,length(tpKeep));
colKeep = colKeep(:);

%check if tracks are in sparse form
sparseForm = issparse(tracks0(1).tracksFeatIndxCG);

%% Sub-sampling

%get number of tracks
numTracks = length(tracks0);

%go over all tracks and sub-sample
tracksSub = tracks0;
if isfield(tracksSub,'aggregState')
    tracksSub = rmfield(tracksSub,'aggregState');
end
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
    
    %find segments that survive the subsampling and those that do not
    indxStay = find(~isnan(segSEL(:,3)));
    indxGone = setdiff((1:numSeg)',indxStay);
    if isempty(indxGone)
        indxGone  = [];
    end
    
    %get the track's sequence of events
    seqOfEvents = tracks0(iTrack).seqOfEvents;
    
    %assign new start and end times - merges will need an addition of one,
    %done two steps down
    for iSeg = 1 : numSeg
        rowsSeg = find(seqOfEvents(:,3)==iSeg);
        seqOfEvents(rowsSeg(1),1) = segSEL(iSeg,1);
        seqOfEvents(rowsSeg(2),1) = segSEL(iSeg,2);
    end
    
    %replace segments that did not survive the subsampling with NaN in both
    %the 3rd and 4th column - this means that some splits followed by a
    %merge might not survive the subsampling
    for iSeg = indxGone'
        seqOfEvents(seqOfEvents(:,3)==iSeg,3) = NaN;
        seqOfEvents(seqOfEvents(:,4)==iSeg,4) = NaN;
    end
    
    %for surviving merges, add one to end time to follow convention
    rowsSeg = find(seqOfEvents(:,2)==2 & ~isnan(seqOfEvents(:,4)));
    seqOfEvents(rowsSeg,1) = seqOfEvents(rowsSeg,1) + 1;
    
    %find merges and splits that involve the same track segment and happen at the exact same time
    %these result from associations that last shorter than the sub-sampling time
    timeSegmentSplit = seqOfEvents(seqOfEvents(:,2)==1&~isnan(seqOfEvents(:,4)),[1 4]);
    timeSegmentMerge = seqOfEvents(seqOfEvents(:,2)==2&~isnan(seqOfEvents(:,4)),[1 4]);
    timeSegmentCommon = intersect(timeSegmentSplit,timeSegmentMerge,'rows');
    
    %go over these merges and splits and cancel them, as basically the
    %merge and split cancel each other out
    %add splitting segment to list of segments not surviving
    for iMS = 1 : size(timeSegmentCommon,1)
        
        %get the merging and splitting segments to be stitched back
        %together
        timeMS = timeSegmentCommon(iMS,1);
        segMS  = timeSegmentCommon(iMS,2);
        segMerge = seqOfEvents( (seqOfEvents(:,1)==timeMS & seqOfEvents(:,4)==segMS & seqOfEvents(:,2)==2) , 3 );
        segSplit = seqOfEvents( (seqOfEvents(:,1)==timeMS & seqOfEvents(:,4)==segMS & seqOfEvents(:,2)==1) , 3 );
        
        %stitch them back together and cancel the merge and split
        trackedFeatureIndx(segMerge,timeMS:end) = trackedFeatureIndx(segSplit,timeMS:end);
        trackedFeatureIndx(segSplit,:) = 0;
        trackedFeatureInfo(segMerge,8*(timeMS-1)+1:end) = trackedFeatureInfo(segSplit,8*(timeMS-1)+1:end);
        trackedFeatureInfo(segSplit,8*(timeMS-1)+1:end) = NaN;
        seqOfEvents(seqOfEvents(:,1)==timeMS & seqOfEvents(:,4)==segMS,3:4) = NaN;
        
        %update segment numbers in sequence of events
        seqOfEvents(seqOfEvents(:,3)==segSplit,3) = segMerge;
        seqOfEvents(seqOfEvents(:,4)==segSplit,4) = segMerge;
        
        %add splitting segment to list of those not surviving
        indxGone = [indxGone; segSplit]; %#ok<AGROW>
        indxStay = setdiff(indxStay,segSplit);
        
    end
    
    %if a first segment merges with a second segment, which itself merges
    %with a third segment at the same time, replace the second segment with
    %the third
    timeSegmentEndMerge = seqOfEvents(seqOfEvents(:,2)==2&~isnan(seqOfEvents(:,4)),[1 3 4]);
    timeSegmentCommon = intersect(timeSegmentEndMerge(:,[1 3]),timeSegmentEndMerge(:,[1 2]),'rows');
    for iM = 1 : size(timeSegmentCommon,1)
        segmentEnd = timeSegmentCommon(iM,2);
        segmentReplace = seqOfEvents(seqOfEvents(:,2)==2&seqOfEvents(:,3)==segmentEnd,4);
        seqOfEvents(seqOfEvents(:,2)==2&seqOfEvents(:,4)==segmentEnd,4) = segmentReplace;
    end
    
    %if a first segment splits from a second segment, which itself splits
    %from a third segment at the same time, replace the second segment with
    %the third
    timeSegmentStartSplit = seqOfEvents(seqOfEvents(:,2)==1&~isnan(seqOfEvents(:,4)),[1 3 4]);
    timeSegmentCommon = intersect(timeSegmentStartSplit(:,[1 3]),timeSegmentStartSplit(:,[1 2]),'rows');
    for iS = 1 : size(timeSegmentCommon,1)
        segmentStart = timeSegmentCommon(iS,2);
        segmentReplace = seqOfEvents(seqOfEvents(:,2)==1&seqOfEvents(:,3)==segmentStart,4);
        seqOfEvents(seqOfEvents(:,2)==1&seqOfEvents(:,4)==segmentStart,4) = segmentReplace;
    end
    
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
    
    %store the sub-sampled compound track
    tracksSub(iTrack).tracksFeatIndxCG = trackedFeatureIndx;
    tracksSub(iTrack).tracksCoordAmpCG = trackedFeatureInfo;
    tracksSub(iTrack).seqOfEvents = seqOfEvents;
    
end

%% ~~~ the end ~~~

