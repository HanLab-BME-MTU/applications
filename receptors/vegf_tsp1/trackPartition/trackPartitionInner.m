function frameIntersection = trackPartitionInner(tracksCoord,seq,maskSub,minTrackLength)

% Number of tracks in this compound track
nTracks = size(tracksCoord,1); 

lengthCompoundTrack = size(tracksCoord,2)/8;

% Find sequence start (start of earliest track)
seqStart = min(seq((seq(:,2) == 1),1));
% Find sequence end (end of latest track)
seqEnd = max(seq((seq(:,2) == 2),1));
tSeq = (seqStart:seqEnd)';

% Store result here
frameIntersection = ones(nTracks,lengthCompoundTrack);



% Iterate through the tracks inside this compound track
for j = 1:nTracks  
    track = tracksCoord(j,:);
    intersectionTemp = frameIntersection(j,:);

    % Find track start and end time from seqOfEvents
    seqTrackStart = (seq(:,2) == 1)&(seq(:,3) == j);
    seqTrackEnd = (seq(:,2) == 2)&(seq(:,3) == j);
    trackStart = seq(seqTrackStart,1);
    trackEnd = seq(seqTrackEnd,1);
    indStart = find(tSeq == trackStart);
    indStop = find(tSeq == trackEnd);

    % Get list of positions (in x and y coords) from the track info
    x = track(1:8:end)';
    y = track(2:8:end)';


    % Interpolate positions during gaps
    if sum(isnan(x)) > 0
        [x(indStart:indStop),y(indStart:indStop)] ...
            = gapInterp(x(indStart:indStop),y(indStart:indStop));
    end

    % Now we have a list of positions at each time:
    % x(t1) y(t1) t1
    % x(t2) y(t2) t2
    % ...
    trackPos = [x,y,tSeq];

    % Find rows in this position list that also appear in the masked
    % coordinate list
    [~,iIntersect,~] = intersect(round(trackPos),maskSub,'rows');
    intersectionTemp(iIntersect) = 100;
    % Pad the beginning and end with zeros where the track does not
    % exist
    if indStart ~= 1
        intersectionTemp(1:indStart) = 0;
    end
    if indStop ~= lengthCompoundTrack
        intersectionTemp(indStop:end) = 0;
    end

    % Merge short tracks
    % should enforce minTrackLength being an odd integer
    if minTrackLength > 0 
        intersectionTemp = mergeShortTracks(intersectionTemp,minTrackLength);
    end

    frameIntersection(j,:) = intersectionTemp;
end
end