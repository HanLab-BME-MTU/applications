function tracksPart = cutTracks(tracks,frameIntersectionList)
%CUTTRACKS Using the intersection information from trackPartitionInner, cut
%up compound tracks each time they enter or exit colocalization

nCompoundTracks = size(tracks,1);

% Store results here
tracksPart = struct('tracksFeatIndxCG',[],'tracksCoordAmpCG',[],...
    'seqOfEvents',[],'isInside',[],'originCompoundTrack',[]);
iRow = 1;

for iCompoundTrack = 1:nCompoundTracks
    % Get existing track data
    tracksIndOld = tracks(iCompoundTrack).tracksFeatIndxCG;
    tracksCoordOld = tracks(iCompoundTrack).tracksCoordAmpCG;
    seqOld = tracks(iCompoundTrack).seqOfEvents;
    
    % The frame the track starts on (used later to convert between relative
    % frame number (frame index within the track) and absolute frame number
    % (frame number within the whole movie))
    frameOffset = min(seqOld(:,1))-1;
    
    frameIntersection = frameIntersectionList{iCompoundTrack};   
    
    sumVector = sum(frameIntersection,1); % Sum down through each track
    lengthCompoundTrack = numel(sumVector);

    segmentsLeft = true;
    % Find the first nonzero value
    segmentStartRel = find(sumVector,1,'first'); 
    segmentStartAbs = segmentStartRel+frameOffset;
    if isempty(segmentStartRel) || (segmentStartRel > lengthCompoundTrack)
        % No nonzero values, so don't cut anything and move on
        segmentsLeft = false;
    end
    while segmentsLeft
        startFrameSum = sumVector(segmentStartRel);
        startFrame = frameIntersection(:,segmentStartRel);
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
        
        % Note: here, the end frame of a segment is the next frame in which
        % it has a different value, i.e. in [1 1 1 100 100], the first
        % segment ends on the fourth frame, not the third. 
                
        % Initialize new sequence of events
        seqNew = [];
        if startFrameSum == 0
            % No tracks exist
            segmentEndRel = find((sumVector(segmentStartRel:end) > 0),1,'first')+segmentStartRel-1;
            if isempty(segmentEndRel)
                % Reached the end of the compound track
                segmentEndRel = lengthCompoundTrack;
                segmentsLeft = false;
            end
        elseif startFrameSum < 100
            % Track(s) exist, are all outside
            
            % Find where this status changes (at least one track becomes
            % inside (sumVector becomes > 100) or the tracks all end 
            % sumVector == 0))
            segmentEndRel = find((sumVector(segmentStartRel:end) > 100)|... 
                (sumVector(segmentStartRel:end) == 0),1,'first')+segmentStartRel-1;
            
            if isempty(segmentEndRel)
                % Reached the end of the compound track
                segmentEndRel = lengthCompoundTrack+1;
                segmentsLeft = false;
            end
            segmentStartAbs = segmentStartRel+frameOffset;
            segmentEndAbs = segmentEndRel+frameOffset;

            seqNew = newSeqOfEvents(frameIntersection,seqOld,segmentStartRel,...
                segmentStartAbs,segmentEndRel,segmentEndAbs,frameOffset);
            
            % Cut this segment from the track coordinates
            if segmentsLeft == false
                % If reached the end of the compound track, just take the
                % segment all the way until the end
                tracksIndNew = tracksIndOld(:,segmentStartRel:end);
                tracksCoordNew = tracksCoordOld(:,(segmentStartRel-1)*8+1:end);
            else
                tracksIndNew = tracksIndOld(:,segmentStartRel:segmentEndRel-1);
                tracksCoordNew = tracksCoordOld(:,(segmentStartRel-1)*8+1:(segmentEndRel-1)*8);
            end
            
            % Store the segment as a new compound track
            tracksPart(iRow).tracksFeatIndxCG = tracksIndNew;
            tracksPart(iRow).tracksCoordAmpCG = tracksCoordNew;
            tracksPart(iRow).seqOfEvents = sortrows(seqNew,1);
            tracksPart(iRow).isInside = false;
            tracksPart(iRow).originCompoundTrack = iCompoundTrack;
            iRow = iRow+1;
        elseif mod(startFrameSum,100) == 0
            % Track(s) exist, all inside
            
            % Find where this status changes (at least one track becomes
            % outside (sumVector becomes not a multiple of 100) or the 
            % tracks all end (sumVector == 0))
            segmentEndRel = find((mod(sumVector(segmentStartRel:end),100) ~= 0)|... 
                (sumVector(segmentStartRel:end) == 0),1,'first')+segmentStartRel-1;

            if isempty(segmentEndRel)
                % Reached the end of the compound track
                segmentEndRel = lengthCompoundTrack+1;
                segmentsLeft = false;
            end
            segmentStartAbs = segmentStartRel+frameOffset;
            segmentEndAbs = segmentEndRel+frameOffset;
            

            seqNew = newSeqOfEvents(frameIntersection,seqOld,segmentStartRel,...
                segmentStartAbs,segmentEndRel,segmentEndAbs,frameOffset);
            
            % Cut this segment from the track coordinates
            if segmentsLeft == false
                % If reached the end of the compound track, just take the
                % segment all the way until the end
                tracksIndNew = tracksIndOld(:,segmentStartRel:end);
                tracksCoordNew = tracksCoordOld(:,(segmentStartRel-1)*8+1:end);
            else
                tracksIndNew = tracksIndOld(:,segmentStartRel:segmentEndRel-1);
                tracksCoordNew = tracksCoordOld(:,(segmentStartRel-1)*8+1:(segmentEndRel-1)*8);
            end
            
            % Store the segment as a new compound track
            tracksPart(iRow).tracksFeatIndxCG = tracksIndNew;
            tracksPart(iRow).tracksCoordAmpCG = tracksCoordNew;
            tracksPart(iRow).seqOfEvents = sortrows(seqNew,1);
            tracksPart(iRow).isInside = true;
            tracksPart(iRow).originCompoundTrack = iCompoundTrack;
            iRow = iRow+1;
        else
            % Some tracks inside and some tracks outside
            
            % Find where this status changes (find the next column of
            % frameIntersection that is different from the segment start
            % column
            segmentEndRel = [];
            for c = (segmentStartRel+1):size(frameIntersection,2)
                column = frameIntersection(:,c);
                if ~isequal(column(:),startFrame(:))
                    segmentEndRel = c;
                    break
                end
            end
            
            if isempty(segmentEndRel)
                % Reached the end of the compound track
                segmentEndRel = lengthCompoundTrack+1;
                segmentsLeft = false;
            end
            segmentStartAbs = segmentStartRel+frameOffset;
            segmentEndAbs = segmentEndRel+frameOffset;
            
            % Deal with inside tracks
            seqNew = newSeqOfEvents(frameIntersection == 100,seqOld,segmentStartRel,...
                segmentStartAbs,segmentEndRel,segmentEndAbs,frameOffset);
            
            
            % Cut this segment from the track coordinates
            if segmentsLeft == false
                % If reached the end of the compound track, just take the
                % segment all the way until the end
                tracksIndNew = tracksIndOld(:,segmentStartRel:end);
                tracksCoordNew = tracksCoordOld(:,(segmentStartRel-1)*8+1:end);
            else
                tracksIndNew = tracksIndOld(:,segmentStartRel:segmentEndRel-1);
                tracksCoordNew = tracksCoordOld(:,(segmentStartRel-1)*8+1:(segmentEndRel-1)*8);
            end
            
            % But keep only the inside tracks
            insideTracksInd = sum(frameIntersection(:,segmentStartRel:segmentEndRel-1)==100,2)' > 0;
            outsideTracksInd = ~insideTracksInd;
            tracksIndNew(outsideTracksInd,:) = 0;
            tracksCoordNew(outsideTracksInd,:) = NaN;
            
            % Store the segment as a new compound track
            tracksPart(iRow).tracksFeatIndxCG = tracksIndNew;
            tracksPart(iRow).tracksCoordAmpCG = tracksCoordNew;
            tracksPart(iRow).seqOfEvents = sortrows(seqNew,1);
            tracksPart(iRow).isInside = true;
            tracksPart(iRow).originCompoundTrack = iCompoundTrack;
            iRow = iRow+1;
            
            % Deal with outside tracks
            seqNew = newSeqOfEvents(frameIntersection == 1,seqOld,segmentStartRel,...
                segmentStartAbs,segmentEndRel,segmentEndAbs,frameOffset);
            
            % Cut this segment from the track coordinates
            if segmentsLeft == false
                % If reached the end of the compound track, just take the
                % segment all the way until the end
                tracksIndNew = tracksIndOld(:,segmentStartRel:end);
                tracksCoordNew = tracksCoordOld(:,(segmentStartRel-1)*8+1:end);
            else
                tracksIndNew = tracksIndOld(:,segmentStartRel:segmentEndRel-1);
                tracksCoordNew = tracksCoordOld(:,(segmentStartRel-1)*8+1:(segmentEndRel-1)*8);
            end
            
            % But keep only the outside tracks
            tracksIndNew(insideTracksInd,:) = 0;
            tracksCoordNew(insideTracksInd,:) = NaN;
            
            % Store the segment as a new compound track
            tracksPart(iRow).tracksFeatIndxCG = tracksIndNew;
            tracksPart(iRow).tracksCoordAmpCG = tracksCoordNew;
            tracksPart(iRow).seqOfEvents = sortrows(seqNew,1);
            tracksPart(iRow).isInside = false;
            tracksPart(iRow).originCompoundTrack = iCompoundTrack;
            iRow = iRow+1;
        end
        segmentStartRel = segmentEndRel;
    end
end
tracksPart = tracksPart';
end

function seqNew = newSeqOfEvents(frameIntersection,seqOld,segmentStartRel,...
        segmentStartAbs,segmentEndRel,segmentEndAbs,frameOffset)
% Find which tracks exist during this segment by summing their
% frameIntersection values up and seeing which sums are >0
% (non-existant tracks are all 0's)
if (segmentEndRel-segmentStartRel) > 0
    validTracksInd = find(sum(frameIntersection(:,segmentStartRel:segmentEndRel-1),2))';
else 
    % segment is only one frame
    validTracksInd = find(frameIntersection(:,segmentStartRel))';
end
nValidTracks = numel(validTracksInd);

seqNew = [];
% Create new seqOfEvents
for iTrack = 1:nValidTracks

    
    % Original track number in seqOfEvents
    trackNumber = validTracksInd(iTrack);

    % Frame at which this track starts or ends
    if (segmentEndRel-segmentStartRel) > 0
        trackStartRel = find(frameIntersection(trackNumber,segmentStartRel:segmentEndRel-1),1,'first');
        trackStartAbs = trackStartRel+segmentStartAbs-1;
        trackEndRel = find(frameIntersection(trackNumber,trackStartRel+segmentStartRel:segmentEndRel-1) == 0,1,'first')+trackStartRel;
        if isempty(trackEndRel)
            trackEndAbs = segmentEndAbs;
        else
            trackEndAbs = trackEndRel+segmentStartAbs-1;
        end
    else
        % segment is only one frame
        trackStartAbs = segmentStartAbs;
        trackEndAbs = segmentEndAbs;
    end

    % Find existing start and end events for this track
    eventStartInd = (seqOld(:,1) == trackStartAbs) & ...
        (seqOld(:,2) == 1) & (seqOld(:,3) == trackNumber);
    if sum(eventStartInd) > 0
        % Use the existing event
        eventStart = seqOld(eventStartInd,:);
    else
        % Make a new track start event
        eventStart = [trackStartAbs,1,trackNumber,NaN];
    end

    eventEndInd = (seqOld(:,1) == trackEndAbs) & ...
        (seqOld(:,2) == 2) & (seqOld(:,3) == trackNumber);
    if sum(eventEndInd) > 0
        % Use the existing event
        eventEnd = seqOld(eventEndInd,:);
    else
        % Make a new track end event
        eventEnd = [trackEndAbs,2,trackNumber,NaN];
    end
    % If the track ends at the end of the segment, make the
    % end frame number the last frame of the segment
    if trackEndAbs == segmentEndAbs;
        eventEnd(1) = segmentEndAbs-1;
    end
    seqNew = [seqNew;eventStart;eventEnd];
end
end
