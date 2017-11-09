function frameIntersection = trackPartitionInner(xCoord,yCoord,startFrame,...
    endFrame,segStartFrame,segEndFrame,maskSub,minTrackLength)
%TRACKPARTITIONINNER Does the meat of the computation for trackPartition
%
%   frameIntersection = trackPartitionInner(xCoord,yCoord,startFrame,...
%            endFrame,segStartFrame,segEndFrame,maskSub,minTrackLength)
%
%   Inputs:
%       xCoord,yCoord:      nSegments x nFrames array of x and y coordinates 
%                           from the compound track
%
%       startFrame:         frame at which the compound track starts
%
%       endFrame:           frame at which the compound track ends
%
%       segStartFrame:      vector of length nSegments containing the start
%                           frame for each segment
%
%       segEndFrame:        vector of length nSegments containing the end
%                           frame for each segment
%
%       maskSub:            3 column array of masked location subscript
%                           indices, yielded by performing find(mask)
%       
%       minTrackLength:     minimum number of frames particle must spend
%                           outside or inside colocalization before its
%                           track is cut into a new segment. Must be an odd
%                           integer.
%
%   Output:
%       frameIntersection:  nSegments x nFrames array with the following
%                           possible values:
%                           0 - segment does not exist
%                           1 - segment is "outside"
%                           100 - segment is "inside"
%
%Kevin Nguyen, July 2016

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
    
    % Find the index in the compound track at which the current segment
    % starts
    segStartInd = find(times == segStartFrame(j));
    segEndInd = find(times == segEndFrame(j));

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
    intersectionTemp(iIntersect) = 100; % Mark intersections with value 100
    
    % Pad the beginning and end with zeros where the track does not
    % exist
    if segStartInd ~= 1
        intersectionTemp(1:segStartInd-1) = 0;
    end
    if segEndInd ~= lengthTrack
        intersectionTemp(segEndInd+1:end) = 0;
    end

    % Merge short tracks
    if minTrackLength > 0 
        intersectionTemp = mergeShortTracks(intersectionTemp,minTrackLength);
    end

    frameIntersection(j,:) = intersectionTemp;
end
end