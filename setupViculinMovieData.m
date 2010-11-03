function movieData = setupViculinMovieData(path,params)

% Analysis directory
movieData.analysisDirectory = fullfile(path, 'ch488', 'analysis');

% Image directory
movieData.imageDirectory = path;
movieData.channelDirectory = cellfun(@(x) fullfile(x, 'roi'), ...
    params.channelDirectory, 'UniformOutput', false);

% Get the number of images
nImages = cellfun(@(channelPath) numel(dir([movieData.imageDirectory ...
    filesep channelPath filesep '*.tif'])), movieData.channelDirectory);
assert(all(nImages(:) == nImages(1)));
movieData.nImages = nImages(1);

% Load physical parameter from
filename = fullfile(movieData.imageDirectory, 'ch560', 'analysis', ...
    'fsmPhysiParam.mat');

if exist(filename, 'file')
    load(filename);
    
    if fsmPhysiParam.pixelSize ~= params.pixelSize
        errMsg = MException(['pixel size in qFSM differs from value in '...
            'batch params (SKIPPING).']);
        throw(errMsg);
    end
    
    if fsmPhysiParam.frameInterval ~= params.timeInterval
        errMsg = MException(['time interval in qFSM differs from value '...
            'in batch params (SKIPPING).']);
        throw(errMsg);
    end
    
    clear fsmPhysiParam;
end

movieData.pixelSize_nm = params.pixelSize;
movieData.timeInterval_s = params.timeInterval;

% Get the mask directory
movieData.masks.channelDirectory = {''};
movieData.masks.directory = fullfile(movieData.imageDirectory, 'ch560', ...
    'analysis', 'edge', 'cell_mask');
if exist(movieData.masks.directory, 'dir')
    movieData.masks.n = numel(dir([movieData.masks.directory filesep '*.tif']));
    movieData.masks.status = 1;
else
    movieData.masks.status = 0;
end
