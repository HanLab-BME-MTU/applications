function movieData = getMovieParticleDetection(movieData, iChannel, sigmaPSF, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.KeepUnmatched = true;
ip.addRequired('movieData');
ip.addRequired('iChannel', @isscalar);
ip.addRequired('sigmaPSF', @isscalar);
ip.addParamValue('batchMode', true, @islogical);

ip.parse(movieData, iChannel, sigmaPSF, varargin{:});
batchMode = ip.Results.batchMode;
unmatched = ip.Unmatched;
unmatched = cellfun(@(field) {field, unmatched.(field)}, ...
    fieldnames(unmatched), 'UniformOutput', false);
unmatched = horzcat(unmatched{:});

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
    'sigmaX',[],...
    'sigmaY',[],...
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
    
    featuresInfo(iFrame) = cometDetection(ima, mask, sigmaPSF, unmatched{:});
        
    if ~batchMode && ishandle(h)
        waitbar(iFrame / nFrames, h);
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
