function tracksPart = trackPartition(tracksInput,mask,imageFrameStart,imageFrameEnd)
%TRACKSPARTITION Finds and labels frames of tracks during which a particle
%intersects (colocates) with a masked area 

imSize = size(mask); % [rowSize,colSize,tSize]
nFrames = imageFrameEnd-imageFrameStart+1;

% Subscript index list for all masked areas in all frames
maskSub = cell(nFrames,1);
for f = imageFrameStart:imageFrameEnd
    maskInd = find(mask(:,:,f));
    maskSubFrame = zeros(numel(maskInd),2);
    [maskSubFrame(:,2),maskSubFrame(:,1)] = ind2sub(imSize,maskInd);
    maskSub{f} = [maskSubFrame,f*ones(numel(maskInd),1)]; % [x(:),y(:),t]
end
maskSub = cell2mat(maskSub);

nCompoundTracks = size(tracksInput,1);

% Store all results here
partitionEventsList = cell(nCompoundTracks,1);
frameIntersectionList = cell(nCompoundTracks,1);

for i = 1:nCompoundTracks
    % Get some fields from the track struct
    seq = tracksInput(i).seqOfEvents;
    tracksCoord = tracksInput(i).tracksCoordAmpCG;
    
    % Number of tracks in this compound track
    nTracks = size(tracksCoord,1); 
    
    lengthCompoundTrack = size(tracksCoord,2)/8;
        
    % Find sequence start (start of earliest track)
    seqStart = min(seq((seq(:,2) == 1),1));
    % Find sequence end (end of latest track)
    seqEnd = max(seq((seq(:,2) == 2),1));
    tSeq = (seqStart:seqEnd)';
    
    % Store result here
    partitionEvents = cell(nTracks,1);
    frameIntersection = ones(nTracks,lengthCompoundTrack);

    
    
    % Iterate through the tracks inside this compound track
    for j = 1:nTracks  
        track = tracksCoord(j,:);
        
        % Find track start and end time from seqOfEvents
        seqTrackStart = (seq(:,2) == 1)&(seq(:,3) == j);
        seqTrackEnd = (seq(:,2) == 2)&(seq(:,3) == j);
        trackStart = seq(seqTrackStart,1);
        trackEnd = seq(seqTrackEnd,1);
        indStart = find(tSeq == trackStart);
        indStop = find(tSeq == trackEnd);
        
        % Get list of positions (in x and y coords) from the track info
%         [x,y] = parseTrack(track,indStart,indStop);
        x = track(1:8:end)';
        y = track(2:8:end)';
%         tTrack = (trackStart:trackEnd)';
        
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
%         % Create a binary vector of the same length as the track
%         intersection = false(numel(x),1);
        % Mark intersection frames
%         intersection(iIntersect) = true;
        frameIntersection(j,iIntersect) = 100;
        % Pad the beginning and end with zeros where the track does not
        % exist
        if indStart ~= 1
            frameIntersection(j,1:indStart) = 0;
        end
        if indStop ~= lengthCompoundTrack
            frameIntersection(j,indStop:end) = 0;
        end
            
%         % Store this binary partition info 
%         frameIntersection{j} = intersection';
        
%         if numel(iIntersect) > 0 % If intersections exist for this track
%             % Are there any intersections left to analyze?
%             intersectionsLeft = true;
%             offset = 0;
%             
%             while intersectionsLeft == true
%                 % Find start and end of next track intersection
%                 intersectStart = find(intersection == 2,1,'first');
%                 intersectEnd = find(intersection(intersectStart:end) == 1,1,'first')+intersectStart-1;
%                 intStartOffset = intersectStart+offset;
%                 intEndOffset = intersectEnd+offset;
%                 
%                 if isempty(intersectEnd)
%                     intersectEnd = indStop;
%                     intEndOffset = indStop;
%                     intersectionsLeft = false;
%                 end
%            
%                 % Store partition info
% %                 newRows = [startInd,3,j,NaN;
% %                             endInd,4,j,NaN];
% % %                 newSeq = insertRow(newSeq,newRows,size(newSeq,1)-2*(nTracks-i)+1);
% %                 newSeq = [newSeq;newRows];
% 
%                 trackEvents = [trackEvents;
%                                  intStartOffset+trackStart-1,1;
%                                  intEndOffset+trackStart-1,2];
%                 % Look for the next intersection starting from the frame
%                 % after this intersection ends
%                 offset = offset+intersectEnd;
%                 intersection = intersection(intersectEnd+1:end);
%                 
%                 % Stop iterating if there are no intersections left
%                 if sum(intersection) == 0
%                     intersectionsLeft = false;
%                 end
%             end 
%             partitionEvents{j} = trackEvents;
%         end
    end
% %     tracksPart(i).seqOfEvents = sortrows(newSeq,1);
%     partitionEventsList{i} = partitionEvents;
    frameIntersectionList{i} = frameIntersection;
end

% Cut up the compound tracks at each partition event, i.e. start a new
% compound track whenever a particle enters or exits a masked area

tracksPart = cutTracks(tracksInput,frameIntersectionList);
% Transpose into an nx1 struct instead of 1xn (to match the input track
% struct)
tracksPart = tracksPart';
end

function [x,y] = parseTrack(trackData,start,stop)
    stop = min(8*stop,size(trackData/8,2));
    x = trackData((start-1)*8+1:8:stop)';
    y = trackData((start-1)*8+2:8:stop)';
end

function newArray = insertRow(oldArray,row,insertAt)
    nRows = size(oldArray,1);
    a1 = 1;
    a2 = insertAt-1;
    if insertAt < nRows
        b1 = insertAt;
        b2 = nRows;
        newArray = [oldArray(a1:a2,:);
                    row;
                    oldArray(b1:b2,:)];
    else
        newArray = [oldArray;row];
    end
end