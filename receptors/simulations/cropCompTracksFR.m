function tracksFR = cropCompTracksFR(tracks0,frameRange)
%cropCompTracksFR crops compound tracks to a certain frame range
%
%SYNOPSIS tracksFR = cropCompTracksFR(tracks0,frameRange)
%
%INPUT  tracks0      : Original tracks, in form of output of
%                      trackCloseGapsKalman.
%       frameRange   : Row vector with minimum and maximum frame to retain.
%
%OUTPUT tracksSub    : Same as input tracks, but only within requested
%                      frame range.
%
%Khuloud Jaqaman, February 2013

%% Input

if nargin < 2
    error('cropCompTracksFR: Incorrect number of input arguments!')
end

%% Cropping

%get number of compount tracks
numTracks = length(tracks0);

%copy input tracks to output tracks0
tracksFR = tracks0;

%go over all compound tracks
for iTrack = 1 : numTracks
    
    %get this track's information
    tracksFeatIndx = tracks0(iTrack).tracksFeatIndxCG;
    tracksCoordAmp = tracks0(iTrack).tracksCoordAmpCG;
    seqOfEvents = tracks0(iTrack).seqOfEvents;
    
    %calculate difference between track start frame and first frame of
    %interest
    frameStartDiff = frameRange(1) - seqOfEvents(1,1) + 1;
    
    %if track starts before first frame of interest ... 
    if frameStartDiff > 1
        
        %remove part before first frame
        tracksFeatIndx = tracksFeatIndx(:,frameStartDiff:end);
        tracksCoordAmp = tracksCoordAmp(:,(frameStartDiff-1)*8+1:end);
        
        %modify the sequence of events
        
        %starts - shift them in time
        indx = find( seqOfEvents(:,1) < frameRange(1) & seqOfEvents(:,2) == 1 & isnan(seqOfEvents(:,4)) );
        seqOfEvents(indx,1) = frameRange(1); %#ok<FNDSB>
        
        %splits - convert them to starts and shift them in time
        indx = find( seqOfEvents(:,1) < frameRange(1) & seqOfEvents(:,2) == 1 & ~isnan(seqOfEvents(:,4)) );
        seqOfEvents(indx,1) = frameRange(1);
        seqOfEvents(indx,4) = NaN;
        
        %ends and merges - remove segment completely
        indx = find( (seqOfEvents(:,1) < frameRange(1) & seqOfEvents(:,2) == 2) | ...
            (seqOfEvents(:,1) <= frameRange(1) & seqOfEvents(:,2) == 2 & ~isnan(seqOfEvents(:,4))) );
        segIndx = seqOfEvents(indx,3);
        seqOfEvents(indx,:) = [];
        for iSeg = 1 : length(segIndx)
            seqOfEvents(seqOfEvents(:,3)==segIndx(iSeg),:) = [];
        end
        
    end
    
    %calculate difference between track end frame and last frame of
    %interest
    frameEndDiff = seqOfEvents(end,1) - frameRange(2) + 1;
    
    %if track ends after last frame of itnerest ...
    if frameEndDiff > 1
        
        %remove part after last frame
        numFramesTmp = frameRange(2) - seqOfEvents(1,1) + 1;
        tracksFeatIndx = tracksFeatIndx(:,1:numFramesTmp);
        tracksCoordAmp = tracksCoordAmp(:,1:numFramesTmp*8);
        
        %modify sequence of events
        
        %ends - shift them in time
        indx = find( seqOfEvents(:,1) > frameRange(2) & seqOfEvents(:,2) == 2 & isnan(seqOfEvents(:,4)) );
        seqOfEvents(indx,1) = frameRange(2); %#ok<FNDSB>
        
        %merges - convert them to ends and shift them in time
        indx = find( seqOfEvents(:,1) > frameRange(2) & seqOfEvents(:,2) == 2 & ~isnan(seqOfEvents(:,4)) );
        seqOfEvents(indx,1) = frameRange(2);
        seqOfEvents(indx,4) = NaN;
        
        %starts and splits - remove segment completely
        indx = find( seqOfEvents(:,1) > frameRange(2) & seqOfEvents(:,2) == 1 );
        segIndx = seqOfEvents(indx,3);
        seqOfEvents(indx,:) = [];
        for iSeg = 1 : length(segIndx)
            seqOfEvents(seqOfEvents(:,3)==segIndx(iSeg),:) = [];
        end
        
    end
    
    %put track information back in structure array
    tracksFR(iTrack).tracksFeatIndxCG = tracksFeatIndx;
    tracksFR(iTrack).tracksCoordAmpCG = tracksCoordAmp;
    tracksFR(iTrack).seqOfEvents = seqOfEvents;
        
end



