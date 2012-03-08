function movieData = getMovieParticleDetection(varargin)

% Input must be:
% varargin{1}: the movie data
movieData = varargin{1};
% varargin{2}: the channel index where the detection needs to be performed
iChannel = varargin{2};
% varargin{3}: the detection function handler
detectFunc = varargin{3};
% varargin{4:end-1}: argument of the detection function, except the ima and
% the mask
% varargin{end}: batch mode
batchMode = varargin{end};

%Indicate that particleDetection was started
movieData.particleDetection.status = 0;

movieData.particleDetection.directory = fullfile(movieData.analysisDirectory, 'particleDetection');

% Create output directory
if ~exist(movieData.particleDetection.directory, 'dir')
    mkdir(movieData.particleDetection.directory);
end

movieData.particleDetection.filename = 'featuresInfo.mat';

imagePath = fullfile(movieData.imageDirectory, movieData.channelDirectory{iChannel});
imageFiles = dir([imagePath filesep '*.tif']);

hasMask = checkMovieMasks(movieData);

if hasMask
    maskPaths = cellfun(@(x) fullfile(movieData.masks.directory, x), ...
        movieData.masks.channelDirectory, 'UniformOutput', false);
    maskFiles = cellfun(@(x) dir([x filesep '*.tif']), maskPaths, 'UniformOutput', false);
end

nFrames = numel(imageFiles);

% Get the image dimension
ima = imread(fullfile(imagePath,imageFiles(1).name));
movieData.imSize = size(ima);
clear ima;

if ~batchMode
    h = waitbar(0,'Please wait, particles detection...');
end

featuresInfo(1:nFrames) = struct(...
    'xCoord',[],...
    'yCoord',[],...
    'amp',[],...
    'stdAlong',[],...
    'stdAside',[],...
    'theta',[],...
    'bkg',[]);

for iFrame = 1:nFrames
    % Read image
    ima = imread(fullfile(imagePath, imageFiles(iFrame).name));
    
    % if masks is available, read the mask
    if hasMask
        mask = imread(fullfile(maskPaths{1}, maskFiles{1}(iFrame).name));
    else
        mask = true(size(ima));
    end
    
    featuresInfo(iFrame) = detectFunc(ima, mask, varargin{4:end-1});
        
    if ~batchMode && ishandle(h)
        waitbar(iFrame/nFrames, h)
    end    
end

if ~batchMode && ishandle(h)
    close(h);
end

save([movieData.particleDetection.directory filesep movieData.particleDetection.filename], ...
    'featuresInfo');

movieData.particleDetection.dateTime = datestr(now);
movieData.particleDetection.status = 1;

updateMovieData(movieData);
