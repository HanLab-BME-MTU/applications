function movieData = getMovieSegmentDetection(movieData,iChannel,sigmaPSF,minSize,batchMode)

%Indicate that detection was started
movieData.segmentDetection.status = 0;

%Check that masks has been performed
assert(checkMovieMasks(movieData));

% Read image file list
imagePath = [movieData.imageDirectory filesep movieData.channelDirectory{iChannel}];
imageFiles = dir([imagePath filesep '*.tif']);

% Read mask file list of every channel
nChannels = numel(movieData.masks.channelDirectory);
maskPaths = cellfun(@(channelPath) fullfile(movieData.masks.directory, channelPath), movieData.masks.channelDirectory, 'UniformOutput', false);
maskFiles = cellfun(@(x) dir([x filesep '*.tif']), maskPaths, 'UniformOutput', false);

movieData.segmentDetection.directory = [movieData.analysisDirectory filesep 'segmentDetection'];

if ~exist(movieData.segmentDetection.directory, 'dir')
    mkdir(movieData.segmentDetection.directory);
end

movieData.segmentDetection.filename = 'segmentParams.mat';

nFrames = numel(imageFiles);

segmentParams = cell(nFrames, 1);

if ~batchMode
    h = waitbar(0,'Please wait, focal adhesion detection....');
end

for i = 1:nFrames
    % Read images
    ima = imread(fullfile(imagePath, imageFiles(i).name));
    
    % Read and merge every channel mask
    mask = zeros(size(ima));
    for iChannel = 1:nChannels
        mask = xor(mask, imread(fullfile(maskPaths{iChannel}, maskFiles{iChannel}(i).name)));
    end
    
    % Get initial segment parameters
    segmentParams{i} = subResSegment2DInit(ima,mask,sigmaPSF,minSize);
    
    if ~batchMode && ishandle(h)
        waitbar(i/nFrames, h)
    end
end

if ~batchMode && ishandle(h)
    close(h);
end

save([movieData.segmentDetection.directory filesep movieData.segmentDetection.filename], ...
    'segmentParams');

movieData.segmentDetection.dateTime = datestr(now);
movieData.segmentDetection.status = 1;

updateMovieData(movieData);