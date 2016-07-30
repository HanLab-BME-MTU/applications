function tracksPart = cutTracks(tracks,frameIntersectionList)
%CUTTRACKS Using the intersection information from trackPartitionInner, cut
%up compound tracks each time they enter or exit colocalization

nCompoundTracks = size(tracks,1);

% Store results here
tracksPart = struct('tracksFeatIndxCG',[],'tracksCoordAmpCG',[],...
    'seqOfEvents',[],'isInside',[],'originCompoundTrack',[]);
iRow = 1;

for iCompoundTrack = 1:nCompoundTracks
    tracksIndOld = tracks(iCompoundTrack).tracksFeatIndxCG;
    tracksCoordOld = tracks(iCompoundTrack).tracksCoordAmpCG;
    seqOld = tracks(iCompoundTrack).seqOfEvents;
    
    % The frame the track starts on
    frameOffset = min(seqOld(:,1))-1;
    
    frameIntersection = frameIntersectionList{iCompoundTrack};   
    
    sumVector = sum(frameIntersection,1); % Sum down through each track
    lengthCompoundTrack = numel(sumVector);

    segmentsLeft = true;
    segmentStart = find(sumVector,1,'first'); % first nonzero value
    if isempty(segmentStart) || (segmentStart > lengthCompoundTrack)
        segmentsLeft = false;
    end
    while segmentsLeft
        startFrameSum = sumVector(segmentStart);
        startFrame = frameIntersection(:,segmentStart);
        % Remember that frameIntersection values can be 
        %   0       track doesn't exist
        %   1       track is outside masked area
        %   100     track is inside masked area
        % Consequently, sumVector values can be 
        %   0       no tracks exist (might happen due to short tracks
        %               getting removed)
        %   n<100   n track(s) exist, are all outside
        %   100*m   m track(s) exist, are all inside
        %   100*m+n m tracks are inside, n tracks are outside
                
        % Create new sequence of events
        seqNew = [];
        if startFrameSum == 0
            % No tracks exist
            segmentEnd = find((sumVector(segmentStart:end) > 0),1,'first')+segmentStart-2;
            if isempty(segmentEnd)
                segmentEnd = lengthCompoundTrack;
                segmentsLeft = false;
            end
        elseif startFrameSum < 100
            % Track(s) exist, all outside (sumVector value is n < 100)
            
            % Find where sumVector value is no longer 0 < n < 100
            segmentEnd = find((sumVector(segmentStart:end) >= 100) | (sumVector(segmentStart:end) == 0),1,'first')+segmentStart-2;
            if isempty(segmentEnd)
                segmentEnd = lengthCompoundTrack;
                segmentsLeft = false;
            end
            validTracks = find(sum(frameIntersection(:,segmentStart:segmentEnd),2))';
            nTracks = numel(validTracks);
            
            tracksIndNew = tracksIndOld(validTracks,segmentStart:segmentEnd);
            tracksCoordNew = tracksCoordOld(validTracks,(segmentStart-1)*8+1:segmentEnd*8);
            
            
            for iTrack = 1:nTracks
                % Find where track starts and ends being non-zero in
                % frameIntersection
                trackStart = find(frameIntersection(validTracks(iTrack),segmentStart:segmentEnd),1,'first')+segmentStart-1+frameOffset;
                trackEnd = find(frameIntersection(validTracks(iTrack),trackStart:segmentEnd)==0,1,'first')+trackStart-1+frameOffset;
                if isempty(trackStart)
                    trackStart = segmentStart+frameOffset;
                end
                if isempty(trackEnd)
                    trackEnd = segmentEnd+frameOffset;
                end
                
                % Find existing events at these start and end times
                indStartEvent = (seqOld(:,1) == trackStart)...
                    &(seqOld(:,2) == 1)&(seqOld(:,3) == validTracks(iTrack));
                indEndEvent = (seqOld(:,1) == trackEnd+1)...
                    &(seqOld(:,2) == 2)&(seqOld(:,3) == validTracks(iTrack));
                if sum(indStartEvent(:)) > 0
                    % There is an existing start event
                    startEvent = seqOld(indStartEvent,:);
                    startEvent(3) = iTrack;
                    if ~isnan(startEvent(4))
                        % Change the split track number
                        try
                            startEvent(4) = findTracks(validTracks == startEvent(4));
                        catch
                            startEvent(4) = 0;
                        end
                    end
                else
                    % No start event, create a new one
                    startEvent = [trackStart,1,iTrack,NaN];
                end
                if sum(indEndEvent(:)) > 0
                    % There is an existing end event
                    endEvent = seqOld(indEndEvent,:);
                    endEvent(3) = iTrack;
                    endEvent(1) = trackEnd;
                    if ~isnan(endEvent(4))
                        % Change the split track number
                        try
                            endEvent(4) = findTracks(validTracks == endEvent(4));
                        catch
                            endEvent(4) = 0;
                        end
                    end
                else
                    % No end event, create a new one
                    endEvent = [trackEnd,2,iTrack,NaN];
                end
                seqNew = [seqNew;startEvent;endEvent];
                isInsideNew = false; 
            end
            
            tracksPart(iRow).tracksFeatIndxCG = tracksIndNew;
            tracksPart(iRow).tracksCoordAmpCG = tracksCoordNew;
            tracksPart(iRow).seqOfEvents = sortrows(seqNew,1);
            tracksPart(iRow).isInside = isInsideNew;
            tracksPart(iRow).originCompoundTrack = iCompoundTrack;

            iRow = iRow+1;
        elseif mod(startFrameSum,100) == 0
            % Track(s) exist, all inside (sumVector value is 100*m)
            
            % Find where frameIntersection value next changes
            for k = 1:nTracks
                frameIntersectionSame = frameIntersection(k,segmentStart:segmentEnd) == startFrame(k);
            end
            segmentEnd = find(sum(frameIntersectionSame,1)==nTracks,1,'first')+segmentStart-2;
            if isempty(segmentEnd)
                segmentEnd = lengthCompoundTrack;
                segmentsLeft = false;
            end
            validTracks = find(sum(frameIntersection(:,segmentStart:segmentEnd),2))';
            nTracks = numel(validTracks);
            
            tracksIndNew = tracksIndOld(validTracks,segmentStart:segmentEnd);
            tracksCoordNew = tracksCoordOld(validTracks,(segmentStart-1)*8+1:segmentEnd*8);
            
            for iTrack = 1:nTracks
                % Find where track starts and ends being non-zero in
                % frameIntersection
                trackStart = find(frameIntersection(validTracks(iTrack),segmentStart:segmentEnd),1,'first')+segmentStart-1+frameOffset;
                trackEnd = find(frameIntersection(validTracks(iTrack),trackStart:segmentEnd)==0,1,'first')+trackStart-1+frameOffset;
                if isempty(trackStart)
                    trackStart = segmentStart+frameOffset;
                end
                if isempty(trackEnd)
                    trackEnd = segmentEnd+frameOffset;
                end
                
                % Find existing events 
                indStartEvent = (seqOld(:,1) == trackStart)...
                    &(seqOld(:,2) == 1)&(seqOld(:,3) == validTracks(iTrack));
                indEndEvent = (seqOld(:,1) == trackEnd+1)...
                    &(seqOld(:,2) == 2)&(seqOld(:,3) == validTracks(iTrack));
                if sum(indStartEvent(:)) > 0
                    % There is an existing start event
                    startEvent = seqOld(indStartEvent,:);
                    startEvent(3) = iTrack;
                    if ~isnan(startEvent(4))
                        % Change the split track number
                        try
                            startEvent(4) = find(validTracks == startEvent(4));
                        catch
                            startEvent(4) = 0;
                        end
                    end
                else
                    % No start event, create a new one
                    startEvent = [trackStart,1,iTrack,NaN];
                end
                if sum(indEndEvent(:)) > 0
                    % There is an existing end event
                    endEvent = seqOld(indEndEvent,:);
                    endEvent(3) = iTrack;
                    endEvent(1) = trackEnd;
                    if ~isnan(endEvent(4))
                        % Change the split track number
                        try
                            endEvent(4) = find(validTracks == endEvent(4));
                        catch
                            endEvent(4) = 0;
                        end
                    end
                else
                    % No end event, create a new one
                    endEvent = [trackEnd,2,iTrack,NaN];
                end
                seqNew = [seqNew;startEvent;endEvent];
                isInsideNew = true; 
            end
            
            tracksPart(iRow).tracksFeatIndxCG = tracksIndNew;
            tracksPart(iRow).tracksCoordAmpCG = tracksCoordNew;
            tracksPart(iRow).seqOfEvents = sortrows(seqNew,1);
            tracksPart(iRow).isInside = isInsideNew;
            tracksPart(iRow).originCompoundTrack = iCompoundTrack;

            iRow = iRow+1;
        elseif mod(startFrameSum,100) > 0
            % Track(s) exist, some inside and some outside (sumVector value
            % is 100*m+n)
            
            % Find where frameIntersection value next changes
            for k = 1:nTracks
                frameIntersectionSame = frameIntersection(k,segmentStart:end) == startFrame(k);
            end
            segmentEnd = find(sum(frameIntersectionSame,1)==nTracks,1,'first')+segmentStart-2;
            
%             segmentEnd = find((frameIntersection(:,segmentStart:end) ~= startFrameI),1,'first')+segmentStart-2;
            if isempty(segmentEnd)
                segmentEnd = lengthCompoundTrack;
                segmentsLeft = false;
            end
            validTracks = find(sum(frameIntersection(:,segmentStart:segmentEnd),2))';
            nTracks = numel(validTracks);
            
            insideTrack = max(frameIntersection(validTracks,segmentStart:segmentEnd),[],2) == 100;
            
            % First, deal with all the outside tracks
            tracksIndNew = tracksIndOld(validTracks(~insideTrack),segmentStart:segmentEnd);
            tracksCoordNew = tracksCoordOld(validTracks(~insideTrack),(segmentStart-1)*8+1:segmentEnd*8);
            seqNewTrack = 1;
            for iTrack = 1:nTracks
                if ~insideTrack(iTrack)
                   % Find where track starts and ends being non-zero in
                    % frameIntersection
                    trackStart = find(frameIntersection(validTracks(iTrack),segmentStart:segmentEnd),1,'first')+segmentStart-1+frameOffset;
                    trackEnd = find(frameIntersection(validTracks(iTrack),trackStart:segmentEnd)==0,1,'first')+trackStart-1+frameOffset;
                    if isempty(trackStart)
                        trackStart = segmentStart+frameOffset;
                    end
                    if isempty(trackEnd)
                        trackEnd = segmentEnd+frameOffset;
                    end

                    % Find existing events 
                    indStartEvent = (seqOld(:,1) == trackStart)...
                        &(seqOld(:,2) == 1)&(seqOld(:,3) == validTracks(iTrack));
                    indEndEvent = (seqOld(:,1) == trackEnd+1)...
                        &(seqOld(:,2) == 2)&(seqOld(:,3) == validTracks(iTrack));
                    if sum(indStartEvent(:)) > 0
                        % There is an existing start event
                        startEvent = seqOld(indStartEvent,:);
                        startEvent(3) = seqNewTrack;
                        if ~isnan(startEvent(4))
                        % Change the split track number
                            try
                                startEvent(4) = find(validTracks == startEvent(4));
                            catch
                                startEvent(4) = 0;
                            end
                        end
                    else
                        % No start event, create a new one
                        startEvent = [trackStart,1,seqNewTrack,NaN];
                    end
                    if sum(indEndEvent(:)) > 0
                        % There is an existing end event
                        endEvent = seqOld(indEndEvent,:);
                        endEvent(3) = seqNewTrack;
                        endEvent(1) = trackEnd;
                        if ~isnan(endEvent(4))
                        % Change the split track number
                            try
                                endEvent(4) = find(validTracks == endEvent(4));
                            catch
                                endEvent(4) = 0;
                            end
                        end
                    else
                        % No end event, create a new one
                        endEvent = [trackEnd,2,seqNewTrack,NaN];
                    end
                    seqNew = [seqNew;startEvent;endEvent];
                    isInsideNew = false;

                    seqNewTrack = seqNewTrack+1;
                end
            end
            tracksPart(iRow).tracksFeatIndxCG = tracksIndNew;
            tracksPart(iRow).tracksCoordAmpCG = tracksCoordNew;
            tracksPart(iRow).seqOfEvents = sortrows(seqNew,1);
            tracksPart(iRow).isInside = isInsideNew;
            tracksPart(iRow).originCompoundTrack = iCompoundTrack;

            iRow = iRow+1;
            
            % Now, deal with all the inside tracks
            seqNew = [];
            seqNewTrack = 1;
            tracksIndNew = tracksIndOld(validTracks(insideTrack),segmentStart:segmentEnd);
            tracksCoordNew = tracksCoordOld(validTracks(insideTrack),(segmentStart-1)*8+1:segmentEnd*8);
            for iTrack = 1:nTracks
                if insideTrack(iTrack)
                    % Find where track starts and ends being non-zero in
                    % frameIntersection
                    trackStart = find(frameIntersection(validTracks(iTrack),segmentStart:segmentEnd),1,'first')+segmentStart-1+frameOffset;
                    trackEnd = find(frameIntersection(validTracks(iTrack),segmentStart:segmentEnd)==0,1,'first')+trackStart-1+frameOffset;
                    if isempty(trackStart)
                        trackStart = segmentStart+frameOffset;
                    end
                    if isempty(trackEnd)
                        trackEnd = segmentEnd+frameOffset;
                    end

                    % Find existing events 
                    indStartEvent = (seqOld(:,1) == trackStart)...
                        &(seqOld(:,2) == 1)&(seqOld(:,3) == validTracks(iTrack));
                    indEndEvent = (seqOld(:,1) == trackEnd+1)...
                        &(seqOld(:,2) == 2)&(seqOld(:,3) == validTracks(iTrack));
                    if sum(indStartEvent(:)) > 0
                        % There is an existing start event
                        startEvent = seqOld(indStartEvent,:);
                        startEvent(3) = seqNewTrack;
                        if ~isnan(startEvent(4))
                        % Change the split track number
                            try
                                startEvent(4) = find(validTracks == startEvent(4));
                            catch
                                startEvent(4) = 0;
                            end
                        end
                    end
                    else
                        % No start event, create a new one
                        startEvent = [trackStart,1,seqNewTrack,NaN];
                    end
                    if sum(indEndEvent(:)) > 0
                        % There is an existing end event
                        endEvent = seqOld(indEndEvent,:);
                        endEvent(3) = seqNewTrack;
                        endEvent(1) = trackEnd;
                        if ~isnan(endEvent(4))
                        % Change the split track number
                            try
                                endEvent(4) = find(validTracks == endEvent(4));
                            catch
                                endEvent(4) = 0;
                            end
                        end
                    
                    else
                        % No end event, create a new one
                        endEvent = [trackEnd,2,seqNewTrack,NaN];
                    end
                    seqNew = [seqNew;startEvent;endEvent];
                    isInsideNew = true;

                    seqNewTrack = seqNewTrack+1;
                end
            end
            
            tracksPart(iRow).tracksFeatIndxCG = tracksIndNew;
            tracksPart(iRow).tracksCoordAmpCG = tracksCoordNew;
            tracksPart(iRow).seqOfEvents = sortrows(seqNew,1);
            tracksPart(iRow).isInside = isInsideNew;
            tracksPart(iRow).originCompoundTrack = iCompoundTrack;

            iRow = iRow+1;
    end
    segmentStart = segmentEnd+1;
end
% Transpose into an nx1 struct instead of 1xn (to match the dimensions of 
% the input track struct)
tracksPart = tracksPart';
end

