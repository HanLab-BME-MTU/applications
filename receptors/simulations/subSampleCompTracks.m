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
%REMARKS NOT FULLY TESTED. DOES NOT FULLY WORK. BUT CHECKING IN SO AS NOT TO LOSE.
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

%% Sub-sampling

%get number of tracks
numTracks = length(tracks0);

%go over all tracks and sub-sample
tracksSub = tracks0;
for iTrack = 1 : numTracks
    
    %convert the current track's information to matrix format
    [trackedFeatureInfo,trackedFeatureIndx] = convStruct2MatIgnoreMS(tracks0(iTrack));
    
    %get each track segment's start, end and life time
    segSEL = getTrackSEL(trackedFeatureInfo);
    numSeg = size(segSEL,1);
    
    %find segments that last for less than the sub-sampling factor
    indxShort = find(segSEL(:,3)<subSampFact);
    
    %go over these segments, and for those that do not overlap with any of
    %the surviving time points, shift them in time to the closest surviving time point
    for iSeg = indxShort'
        tpSeg = segSEL(iSeg,1) : segSEL(iSeg,2);
        if isempty(intersect(tpSeg,tpKeep))
            startTP0 = tpSeg(1);
            endTP0 = tpSeg(end);
            startDiffAbs = abs(startTP0 - tpKeep);
            endDiffAbs = abs(endTP0 - tpKeep);
            if min(startDiffAbs) <= min(endDiffAbs)
                subTP = tpKeep(startDiffAbs==min(startDiffAbs));
                startCol0 = (startTP0-1)*8 + 1;
                subCol = (subTP-1)*8 + 1;
                trackedFeatureIndx(iSeg,subTP) = trackedFeatureIndx(iSeg,startTP0);
                trackedFeatureInfo(iSeg,subCol:subCol+7) = trackedFeatureInfo(iSeg,startCol0:startCol0+7);
            else
                subTP = tpKeep(endDiffAbs==min(endDiffAbs));
                endCol0 = (endTP0-1)*8 + 1;
                subCol = (subTP-1)*8 + 1;
                trackedFeatureIndx(iSeg,subTP) = trackedFeatureIndx(iSeg,endTP0);
                trackedFeatureInfo(iSeg,subCol:subCol+7) = trackedFeatureInfo(iSeg,endCol0:endCol0+7);
            end
        end
    end
    
    %keep only the time points aof interest
    trackedFeatureIndx = trackedFeatureIndx(:,tpKeep);
    trackedFeatureInfo = trackedFeatureInfo(:,colKeep);
    
    %get the new start and end time of each track segment
    segSEL = getTrackSEL(trackedFeatureInfo);
    
    %modify the track's sequence of events
    seqOfEvents = tracks0(iTrack).seqOfEvents;
    for iSeg = 1 : numSeg
        seqOfEvents(seqOfEvents(:,2)==1&seqOfEvents(:,3)==iSeg,1) = segSEL(iSeg,1);
        seqOfEvents(seqOfEvents(:,2)==2&seqOfEvents(:,3)==iSeg,1) = segSEL(iSeg,2);
    end
    rowIndx = find(seqOfEvents(:,2)==2&~isnan(seqOfEvents(:,4)));
    seqOfEvents(rowIndx,1) = seqOfEvents(rowIndx,1) + 1;
    
    %store the sub-sampled compound track
    tracksSub(iTrack).tracksFeatIndxCG = trackedFeatureIndx;
    tracksSub(iTrack).tracksCoordAmpCG = trackedFeatureInfo;
    tracksSub(iTrack).seqOfEvents = seqOfEvents;
    
end


