function movieData = makeFAFigure1(movieData,batchMode)

%Indicate that figure creation was started
movieData.figures.status = 0;

% Check that detection has been performed
assert(checkMovieDetection(movieData) == true);

% Check that tracking has been performed
assert(checkMovieTracking(movieData) == true);

movieData.figures.directory = [movieData.channels(1).analysisDirectory ...
    filesep 'figures'];

if ~exist(movieData.figures.directory, 'dir')
    mkdir(movieData.figures.directory);
end

% Load detections
load(fullfile(movieData.detection.directory, movieData.detection.filename));

nFrames = numel(segmentParams); %#ok<USENS>

% Load tracks
load(fullfile(movieData.tracks.directory, movieData.tracks.filename));

% Remove tracks which either
% - lifetime < 2,
% - starts at first frame
% - ends at last frame

trackSEL = getTrackSEL(tracksFinal); %#ok<NODEF>
ind = trackSEL(:,3) >= 2 | trackSEL(:,1) ~= 1 | trackSEL(:,2) ~= nFrames;
tracksFinal = tracksFinal(ind);
trackSEL = trackSEL(ind);

% Compute translocation (distance between FA position at birth and death)
nTracks = numel(tracksFinal);

transLoc = zeros(nTracks,1);

for iTrack = 1:nTracks
    %get track start and end times
    startTime = trackSEL(iTrack,1);
    endTime   = trackSEL(iTrack,2);
    
    x1 = segmentParams{startTime}(tracksFinal(iTrack).tracksFeatIndxCG(1),1);
    x2 = segmentParams{endTime}(tracksFinal(iTrack).tracksFeatIndxCG(end),1);
    
    y1 = segmentParams{startTime}(tracksFinal(iTrack).tracksFeatIndxCG(1),2);
    y2 = segmentParams{ebdTime}(tracksFinal(iTrack).tracksFeatIndxCG(end),2);
    
    transLoc(iTrack) = sqrt((x1 - x2)^2 + (y1 - y2)^2);
end

movieData.figures.dateTime = datestr(now);
movieData.figures.status = 1;

updateMovieData(movieData);
