function movieData = getMovieDetection(movieData,channelIndex,sigmaPSF,minSize,batchMode)

%Indicate that detection was started
movieData.detection.status = 0;

%Check that masks has been performed
assert(checkMovieMasks(movieData));

% Read image file list
imagePath = [movieData.imageDirectory filesep movieData.channelDirectory{channelIndex}];
imageFiles = dir([imagePath filesep '*.tif']);

% Read mask file list of every channel
nChannels = numel(movieData.masks.channelDirectory);
maskPaths = cellfun(@(channelPath) fullfile(movieData.masks.directory, channelPath), movieData.masks.channelDirectory);
maskFiles = cell(nChannels,1);
for iChannel = 1:nChannels
    maskFiles{iChannel} = dir([maskPaths{iChannel} filesep '*.tif']);
end

movieData.detection.directory = [movieData.analysisDirectory filesep 'detection'];

if ~exist(movieData.detection.directory, 'dir')
    mkdir(movieData.detection.directory);
end

movieData.detection.filename = 'segmentParams.mat';

nFrames = numel(imageFiles);

segmentParams = cell(nFrames, 1);

if ~batchMode
    h = waitbar(0,'Please wait, focal adhesion detection....');
end

for i = 1:nFrames
    % Read images
    ima = imread(fullfile(imagePath, imageFiles(i).name));
    
    % Make sure I is type double
    ima = double(ima);

    % Read and merge every channel mask
    mask = zeros(size(ima));
    for iChannel = 1:nChannels
        mask = mask | imread(fullfile(maskPaths{iChannel}, maskFiles{iChannel}(i).name));
    end
    
    % Get initial segment parameters
    segmentParams{i} = getInitialSegmentParams(ima,mask,sigmaPSF,minSize);
    
    if ~batchMode && ishandle(h)
        waitbar(i/nFrames, h)
    end
end

if ~batchMode && ishandle(h)
    close(h);
end

save([movieData.detection.directory filesep movieData.detection.filename], ...
    'segmentParams');

movieData.detection.dateTime = datestr(now);
movieData.detection.status = 1;

updateMovieData(movieData);