function tracksPart = trackPartition(tracksInput,mask,minTrackLength,frameStart,frameEnd)
%TRACKSPARTITION Performs track partitioning, i.e. dividing particle tracks
%into new segments based on when they colocalize with another particle's
%location
%
%   tracksPart = trackPartition(tracksInput,mask,minTrackLength,imageFrameStart,imageFrameEnd)
%
%   Inputs:
%       tracksInput:        track struct, such as one produced by
%                           TrackingProcess. Should contain the fields:
%                           .tracksFeatIndxCG
%                           .tracksCoordAmpCG
%                           .seqOfEvents
%
%       mask:               binary mask image from trackPartitionInit       
%
%       minTrackLength:     minimum number of frames particle must spend
%                           outside or inside colocalization before its
%                           track is cut into a new segment. Also, tracks
%                           shorter than this length will be removed from
%                           the output track struct. minTrackLength must be
%                           an odd integer.
%       
%       frameStart:         movie frame at which to start analysis
%
%       frameEnd:           movie frame at which to stop analysis
%
%   Output:
%       tracksPart:         new track struct in which the original compound
%                           tracks have been cut (split) whenever the tracked
%                           particle enters or exits colocalization. Struct
%                           fields .tracksFeatIndxCG and .tracksCoordAmpCG
%                           are of the same format as that of the input
%                           tracks. The following fields are new/different:
%                           .seqOfEvents:   has two new event types:
%                                           3 - particle enters
%                                           colocalization
%                                           4 - particle exits
%                                           colocalization
%                           .isInside:      'true' if a track is "inside"
%                                           colocalization and 'false' if 
%                                           it is not
%                           .originCompoundTrack:
%                                           number of the original compound
%                                           track in the input tracks from
%                                           which this new compound track
%                                           was cut
%
%Kevin Nguyen, July 2016

% minTrackLength must be odd
assert(mod(minTrackLength,2) == 1,'Minimum track length must be an odd integer')

imSize = size(mask); % [rowSize,colSize,tSize]
nFrames = frameEnd-frameStart+1;

% Create a list of subscript indices for all masked areas in all frames
maskSub = cell(nFrames,1);
for f = frameStart:frameEnd
    maskInd = find(mask(:,:,f));
    maskSubFrame = zeros(numel(maskInd),2);
    [maskSubFrame(:,2),maskSubFrame(:,1)] = ind2sub(imSize,maskInd);
    maskSub{f} = [maskSubFrame,f*ones(numel(maskInd),1)]; % [x(:),y(:),t]
end
maskSub = int8(cell2mat(maskSub));

% Create a TracksHandle object from the input track struct (thanks, Mark)
tracksObj = TracksHandle(tracksInput);
[xCoords,yCoords,startFrames,endFrames,segStartFrames,segEndFrames] = ...
    parseTracksObj(tracksObj);

fprintf('Partitioning tracks... \n')
frameIntersections = parcellfun_progress(@(xCoords,yCoords,startFrames,...
    endFrames,segStartFrames,segEndFrames) ...
    trackPartitionInner(xCoords,yCoords,startFrames,endFrames,...
    segStartFrames,segEndFrames,maskSub,minTrackLength),...
    xCoords,yCoords,startFrames,endFrames,segStartFrames,segEndFrames,...
    'UniformOutput',false);
% pool = gcp;
% job = batch(@run.runFunction,1,{xCoords,yCoords,startFrames,endFrames,segStartFrames,segEndFrames,maskSub,minTrackLength});
% wait(job,'finished');
% frameIntersections = fetchOutputs(job);
% delete(job);
    


% Cut up the compound tracks at each partition event, i.e. start a new
% compound track whenever a particle enters or exits a masked area

tracksPart = cutTracks(tracksInput,frameIntersections);

% Get rid of short compound tracks
nCompoundTracks = numel(tracksPart);
tracksGood = false(nCompoundTracks,1);
for i = 1:nCompoundTracks
    if size(tracksPart(i).tracksFeatIndxCG,2) >= minTrackLength
        tracksGood(i) = true;
    end
end
tracksPart = tracksPart(tracksGood);
tracksPart = tracksPart';
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
    xCoords{i} = int8(trackObj(i).x);
    yCoords{i} = int8(trackObj(i).y);
    startFrames{i} = trackObj(i).startFrame;
    endFrames{i} = trackObj(i).endFrame;
    segStartFrames{i} = trackObj(i).segmentStartFrame;
    segEndFrames{i} = trackObj(i).segmentEndFrame;
end
end