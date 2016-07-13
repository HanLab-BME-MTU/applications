function frameIntersection = trackPartitionInner(xCoord,yCoord,startFrame,...
    endFrame,segStartFrame,segEndFrame,maskSub,minTrackLength)

% Number of segments in this compound track
nSegments = size(xCoord,1); 

% Length of the compound track
lengthTrack = endFrame-startFrame+1;

% Vector of frame numbers
times = (startFrame:endFrame)';

% Store result here
frameIntersection = ones(nSegments,lengthTrack);

% Iterate through the segments inside the compound track
for j = 1:nSegments  
    intersectionTemp = frameIntersection(j,:);

%     % Find track start and end time from seqOfEvents
%     seqTrackStart = (seq(:,2) == 1)&(seq(:,3) == j);
%     seqTrackEnd = (seq(:,2) == 2)&(seq(:,3) == j);
%     trackStart = seq(seqTrackStart,1);
%     trackEnd = seq(seqTrackEnd,1);
%     indStart = find(times == trackStart);
%     indStop = find(times == trackEnd);
    
    % Find the index in the compound track at which the current segment
    % starts
    segStartInd = find(times == segStartFrame(j));
    segEndInd = find(times == segEndFrame(j));
% 
%     % Get list of positions (in x and y coords) from the track info
%     x = track(1:8:end)';
%     y = track(2:8:end)';
    x = xCoord(j,:)';
    y = yCoord(j,:)';


    % Interpolate positions during gaps
    if sum(isnan(xCoord)) > 0
        [x(segStartInd:segEndInd),y(segStartInd:segEndInd)] ...
            = gapInterp(x(segStartInd:segEndInd),y(segStartInd:segEndInd));
    end

    % Now we have a list of coordinates at each time:
    % x(t1) y(t1) t1
    % x(t2) y(t2) t2
    % ...
    coordList = [x,y,times];

    % Find rows in this position list that also appear in the masked
    % coordinate list
    % There may be some faster versions of this intersect function on
    % MATLAB File Exchange
    [~,iIntersect,~] = intersect(round(coordList),maskSub,'rows');
    intersectionTemp(iIntersect) = 100;
    % Pad the beginning and end with zeros where the track does not
    % exist
    if segStartInd ~= 1
        intersectionTemp(1:segStartInd) = 0;
    end
    if segEndInd ~= lengthTrack
        intersectionTemp(segEndInd:end) = 0;
    end

    % Merge short tracks
    if minTrackLength > 0 
        intersectionTemp = mergeShortTracks(intersectionTemp,minTrackLength);
    end

    frameIntersection(j,:) = intersectionTemp;
end
end