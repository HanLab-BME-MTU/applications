function tracksPart = trackPartition(tracksInput,mask,imageFrameStart,imageFrameEnd,minTrackLength)
%TRACKSPARTITION Finds and labels frames of tracks during which a particle
%intersects (colocates) with a masked area 

% minTrackLength must be odd
assert(mod(minTrackLength,2) == 1,'Minimum track length must be an odd integer')

imSize = size(mask); % [rowSize,colSize,tSize]
nFrames = imageFrameEnd-imageFrameStart+1;

% Create a list of subscript indices for all masked areas in all frames
maskSub = cell(nFrames,1);
for f = imageFrameStart:imageFrameEnd
    maskInd = find(mask(:,:,f));
    maskSubFrame = zeros(numel(maskInd),2);
    [maskSubFrame(:,2),maskSubFrame(:,1)] = ind2sub(imSize,maskInd);
    maskSub{f} = [maskSubFrame,f*ones(numel(maskInd),1)]; % [x(:),y(:),t]
end
maskSub = cell2mat(maskSub);

% % Convert input track struct into cell arrays
% cellTemp = struct2cell(tracksInput);
% tracksCoordCell = cellTemp(2,:)';
% seqCell = cellTemp(3,:)';
% 
% fprintf('Partitioning tracks... \n')
% frameIntersections = parcellfun_progress(@(x,y) trackPartitionInner(x,y,maskSub,minTrackLength),...
%     tracksCoordCell,seqCell,'UniformOutput',false);

% Create a TracksHandle object from the input track struct (thanks, Mark)
tracksObj = TracksHandle(tracksInput);
% Collect various properties into cell arrays
% xCoords = tracksObj.x;
% yCoords = tracksObj.y;
% startFrames = tracksObj.startFrame;
% endFrames = tracksObj.endFrame;
% segStartFrames = tracksObj.segmentStartFrame;
% segEndFrames = tracksObj.segmentEndFrame;
[xCoords,yCoords,startFrames,endFrames,segStartFrames,segEndFrames] = ...
    parseTracksObj(tracksObj);

fprintf('Partitioning tracks... \n')
frameIntersections = parcellfun_progress(@(xCoords,yCoords,startFrames,...
    endFrames,segStartFrames,segEndFrames) ...
    trackPartitionInner(xCoords,yCoords,startFrames,endFrames,...
    segStartFrames,segEndFrames,maskSub,minTrackLength),...
    xCoords,yCoords,startFrames,endFrames,segStartFrames,segEndFrames,...
    'UniformOutput',false);

% Cut up the compound tracks at each partition event, i.e. start a new
% compound track whenever a particle enters or exits a masked area

tracksPart = cutTracks(tracksInput,frameIntersections);
% Transpose into an nx1 struct instead of 1xn (to match the dimensions of 
% the input track struct)
tracksPart = tracksPart';

% Get rid of short compound tracks
nCompoundTracks = numel(tracksPart);
tracksGood = false(nCompoundTracks,1);
for i = 1:nCompoundTracks
    if size(tracksPart(i).tracksFeatIndxCG,2) >= minTrackLength
        tracksGood(i) = true;
    end
end
tracksPart = tracksPart(tracksGood);
end

function [xCoords,yCoords,startFrames,endFrames,segStartFrames,...
    segEndFrames] = parseTracksObj(trackObj)
nTracks = size(trackObj,1);
xCoords = cell(nTracks,1);
yCoords = cell(nTracks,1);
startFrames = cell(nTracks,1);
endFrames = cell(nTracks,1);
segStartFrames = cell(nTracks,1);
segEndFrames = cell(nTracks,1);

for i = 1:nTracks
    xCoords{i} = trackObj(i).x;
    yCoords{i} = trackObj(i).y;
    startFrames{i} = trackObj(i).startFrame;
    endFrames{i} = trackObj(i).endFrame;
    segStartFrames{i} = trackObj(i).segmentStartFrame;
    segEndFrames{i} = trackObj(i).segmentEndFrame;
end
end