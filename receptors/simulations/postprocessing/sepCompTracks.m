function tracksOut = sepCompTracks(tracksIn)
%SEPCOMPTRACKS separates tracks not interacting with each other into distinct compound tracks; such issues occur after sub-sampling or time-cropping
%
%SYNOPSIS tracksOut = sepCompTracks(tracksIn)
%
%INPUT  tracksIn     : Original tracks, in form of output of
%                      trackCloseGapsKalman.
%
%OUTPUT tracksOut    : Same as input tracks, just separated into distinct
%                      tracks as needed.
%
%Khuloud Jaqaman, December 2014

%% Input

if nargin < 1
    error('sepCompTracks: Incorrect number of input arguments!')
end

%% Track separation

%get number of tracks
numTracks = length(tracksIn);

%pre-allocate memory for output for efficiency
tracksOut = [tracksIn;tracksIn];

%get number of segments per track
numSegPerTrack = getNumSegmentsPerTrack(tracksIn);
indxTrack = find(numSegPerTrack>1)';

%go over all tracks that have more than one segment
for iTrack = indxTrack

    %get the track's sequence of events
    seqOfEvents = tracksIn(iTrack).seqOfEvents;
    
    %identify groups of segments that interact with each other
    segGroup = zeros(numSegPerTrack(iTrack));
    segGroup(:,1) = (1 : numSegPerTrack(iTrack))';
    interPairs = seqOfEvents(~isnan(seqOfEvents(:,4)),3:4);
    while ~isempty(interPairs)
        groupSeed = interPairs(1,1);
        indxAdd = find(ismember(interPairs(:,1),groupSeed)|ismember(interPairs(:,2),groupSeed));
        while ~isempty(indxAdd)
            groupAdd = interPairs(indxAdd,:);
            interPairs(indxAdd,:) = [];
            groupSeed = unique([groupSeed; groupAdd(:)]);
            indxAdd = find(ismember(interPairs(:,1),groupSeed)|ismember(interPairs(:,2),groupSeed));
        end
        segGroup(groupSeed(1),1:length(groupSeed)) = groupSeed';
        segGroup(ismember(segGroup(:,1),segGroup(groupSeed(1),2:end)'),1) = 0;
    end
    segGroup(segGroup(:,1)==0,:) = [];
    numGroups = size(segGroup,1);
    
    %if there is more than 1 group of interacting segments, save each group
    %in its own compound track
    if numGroups > 1
        for iGroup = 1 : numGroups
            
            %get this group's segments
            indxSeg = segGroup(iGroup,:);
            indxSeg = indxSeg(indxSeg~=0);
            
            %extract information for this group's segments
            tracksFeatIndxCG = tracksIn(iTrack).tracksFeatIndxCG(indxSeg,:);
            tracksCoordAmpCG = tracksIn(iTrack).tracksCoordAmpCG(indxSeg,:);
            indxEvent = ismember(seqOfEvents(:,3),indxSeg) | ismember(seqOfEvents(:,4),indxSeg);
            seqOfEventsTmp = seqOfEvents(indxEvent,:);
            
            %remove extra entries before local start and after local end
            startShift = seqOfEventsTmp(1,1) - seqOfEvents(1,1);
            endShift   = seqOfEvents(end,1) - seqOfEventsTmp(end,1);
            tracksFeatIndxCG = tracksFeatIndxCG(:,startShift+1:end-endShift);
            tracksCoordAmpCG = tracksCoordAmpCG(:,startShift*8+1:end-8*endShift);
            
            %renumber segments to reflect new compound track
            for iStay = 1 : length(indxSeg)
                iSeg = indxSeg(iStay);
                seqOfEventsTmp(seqOfEventsTmp(:,3)==iSeg,3) = iStay;
                seqOfEventsTmp(seqOfEventsTmp(:,4)==iSeg,4) = iStay;
            end
            
            %save in output structure
            if iGroup == 1 %save first group in place of old compound track
                tracksOut(iTrack).tracksFeatIndxCG = tracksFeatIndxCG;
                tracksOut(iTrack).tracksCoordAmpCG = tracksCoordAmpCG;
                tracksOut(iTrack).seqOfEvents = seqOfEventsTmp;
            else %append tracks structure with additional groups
                numTracks = numTracks + 1;
                tracksOut(numTracks).tracksFeatIndxCG = tracksFeatIndxCG;
                tracksOut(numTracks).tracksCoordAmpCG = tracksCoordAmpCG;
                tracksOut(numTracks).seqOfEvents = seqOfEventsTmp;
            end
            
        end
    end
    
end %(for iTrack = indxTrack)

%keep only entries with valid tracks
tracksOut = tracksOut(1:numTracks);

%% ~~~ the end ~~~

