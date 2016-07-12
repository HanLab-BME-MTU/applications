function tracksPart = cutTracks(tracks,frameIntersectionList)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

nCompoundTracks = size(tracks,1);

% Store results here
tracksPart = struct('tracksFeatIndxCG',[],'tracksCoordAmpCG',[],...
    'seqOfEvents',[],'isInside',[],'originCompoundTrack',[]);
iRow = 1;

for iCompoundTrack = 1:nCompoundTracks
    tracksIndOld = tracks(iCompoundTrack).tracksFeatIndxCG;
    tracksCoordOld = tracks(iCompoundTrack).tracksCoordAmpCG;
    seqOld = tracks(iCompoundTrack).seqOfEvents;
    
    frameIntersection = frameIntersectionList{iCompoundTrack};   
    
    sumVector = sum(frameIntersection,1); % Sum down through each track

    segmentsLeft = true;
    segmentStart = find(sumVector,1,'first'); % first nonzero value
    while segmentsLeft
        startFrame = sumVector(segmentStart);
        % Remember that frameIntersection values can be 
        %   0       track doesn't exist
        %   1       track is outside masked area
        %   100     track is inside masked area
        % Consequently, sumVector values can be 
        %   0       no tracks exist (which shouldn't ever happen)
        %   n<100   n track(s) exist, are all outside
        %   100*m   m track(s) exist, are all inside
        %   100*m+n m tracks are inside, n tracks are outside
        validTracks = find(frameIntersection(:,segmentStart))';
        nTracks = numel(validTracks);
        segmentEnd = find(sumVector(segmentStart:end) ~= startFrame,1,'first')+segmentStart-2;
        if isempty(segmentEnd)
            segmentEnd = numel(sumVector);
            segmentsLeft = false;
        end
        tracksIndNew = tracksIndOld(validTracks,segmentStart:segmentEnd);
        tracksCoordNew = tracksCoordOld(validTracks,(segmentStart-1)*8+1:segmentEnd*8);
        
        % Create new sequence of events
        seqNew = [];
        if startFrame == 0
            % No tracks exist
            segmentsLeft = false;
        elseif startFrame < 100
            % Track(s) exist, all outside        
            for iTrack = 1:numel(validTracks)
                % Find existing events 
                indStartEvent = (seqOld(:,1) == segmentStart)&(seqOld(:,3) == validTracks(iTrack));
                indEndEvent = (seqOld(:,1) == segmentEnd)&(seqOld(:,3) == validTracks(iTrack));
                if sum(indStartEvent(:)) > 0
                    % There is an existing start event
                    startEvent = seqOld(indStartEvent,:);
                    startEvent(3) = iTrack;
                else
                    % No start event, create a new one
                    startEvent = [segmentStart,1,iTrack,NaN];
                end
                if sum(indEndEvent(:)) > 0
                    % There is an existing end event
                    endEvent = seqOld(indEndEvent,:);
                    endEvent(3) = iTrack;
                else
                    % No end event, create a new one
                    endEvent = [segmentEnd,2,iTrack,NaN];
                end
                seqNew = [seqNew;startEvent;endEvent];
                isInsideNew = false; 
            end
            
            tracksPart(iRow).tracksFeatIndxCG = tracksIndNew;
            tracksPart(iRow).tracksCoordAmpCG = tracksCoordNew;
            tracksPart(iRow).seqOfEvents = seqNew;
            tracksPart(iRow).isInside = isInsideNew;
            tracksPart(iRow).originCompoundTrack = iCompoundTrack;

            iRow = iRow+1;
        elseif mod(startFrame,100) == 0
            % Track(s) exist, all inside
            for iTrack = 1:numel(validTracks)
                % Find existing events 
                indStartEvent = (seqOld(:,1) == segmentStart)&(seqOld(:,3) == validTracks(iTrack));
                indEndEvent = (seqOld(:,1) == segmentEnd)&(seqOld(:,3) == validTracks(iTrack));
                if sum(indStartEvent(:)) > 0
                    % There is an existing start event
                    startEvent = seqOld(indStartEvent,:);
                    startEvent(3) = iTrack;
                else
                    % No start event, create a new one
                    startEvent = [segmentStart,1,iTrack,NaN];
                end
                if sum(indEndEvent(:)) > 0
                    % There is an existing end event
                    endEvent = seqOld(indEndEvent,:);
                    endEvent(3) = iTrack;
                else
                    % No end event, create a new one
                    endEvent = [segmentEnd,2,iTrack,NaN];
                end
                seqNew = [seqNew;startEvent;endEvent];
                isInsideNew = true; 
            end
            
            tracksPart(iRow).tracksFeatIndxCG = tracksIndNew;
            tracksPart(iRow).tracksCoordAmpCG = tracksCoordNew;
            tracksPart(iRow).seqOfEvents = seqNew;
            tracksPart(iRow).isInside = isInsideNew;
            tracksPart(iRow).originCompoundTrack = iCompoundTrack;

            iRow = iRow+1;
        elseif mod(startFrame,100) > 0
            % Track(s) exist, some inside and some outside
            % First, deal with all the outside tracks
            for iTrack = 1:numel(validTracks)
                if frameIntersection(iTrack,segmentStart) == 1
                    % Find existing events 
                    indStartEvent = (seqOld(:,1) == segmentStart)&(seqOld(:,3) == validTracks(iTrack));
                    indEndEvent = (seqOld(:,1) == segmentEnd)&(seqOld(:,3) == validTracks(iTrack));
                    if sum(indStartEvent(:)) > 0
                        % There is an existing start event
                        startEvent = seqOld(indStartEvent,:);
                        startEvent(3) = iTrack;
                    else
                        % No start event, create a new one
                        startEvent = [segmentStart,1,iTrack,NaN];
                    end
                    if sum(indEndEvent(:)) > 0
                        % There is an existing end event
                        endEvent = seqOld(indEndEvent,:);;
                        endEvent(3) = iTrack;
                    else
                        % No end event, create a new one
                        endEvent = [segmentEnd,2,iTrack,NaN];
                    end
                    seqNew = [seqNew;startEvent;endEvent];
                    isInsideNew = false;
                    
                    tracksPart(iRow).tracksFeatIndxCG = tracksIndNew;
                    tracksPart(iRow).tracksCoordAmpCG = tracksCoordNew;
                    tracksPart(iRow).seqOfEvents = seqNew;
                    tracksPart(iRow).isInside = isInsideNew;
                    tracksPart(iRow).originCompoundTrack = iCompoundTrack;

                    iRow = iRow+1;
                end
            end
            
            % Now, deal with all the inside tracks
            seqNew = [];
            for iTrack = 1:numel(validTracks)
                if frameIntersection(iTrack,segmentStart) == 100
                    % Find existing events 
                    indStartEvent = (seqOld(:,1) == segmentStart)&(seqOld(:,3) == validTracks(iTrack));
                    indEndEvent = (seqOld(:,1) == segmentEnd)&(seqOld(:,3) == validTracks(iTrack));
                    if sum(indStartEvent(:)) > 0
                        % There is an existing start event
                        startEvent = seqOld(indStartEvent,:);
                        startEvent(3) = iTrack;
                    else
                        % No start event, create a new one
                        startEvent = [segmentStart,1,iTrack,NaN];
                    end
                    if sum(indEndEvent(:)) > 0
                        % There is an existing end event
                        endEvent = seqOld(indEndEvent,:);
                        endEvent(3) = iTrack;
                    else
                        % No end event, create a new one
                        endEvent = [segmentEnd,2,iTrack,NaN];
                    end
                    seqNew = [seqNew;startEvent;endEvent];
                    isInsideNew = true;
                    
                    tracksPart(iRow).tracksFeatIndxCG = tracksIndNew;
                    tracksPart(iRow).tracksCoordAmpCG = tracksCoordNew;
                    tracksPart(iRow).seqOfEvents = seqNew;
                    tracksPart(iRow).isInside = isInsideNew;
                    tracksPart(iRow).originCompoundTrack = iCompoundTrack;

                    iRow = iRow+1;
                end
            end
        end
        segmentStart = segmentEnd+1;
    end
end

