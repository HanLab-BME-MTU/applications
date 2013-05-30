%runTracking(data, varargin) tracks CCPs in the movies passed with the 'data' structure.
% This function generates a list of tracks in 'Tracking/trackedFeatures.mat' for each
% data set.
%
% Inputs   
%                  data : list of movies, using the structure returned by loadConditionData.m
%            {settings} : tracker settings. Default: 
%                         loadTrackSettings('Radius', [3 6], 'MaxGapLength', 2);    
%
% Options ('specifier', value)
%           'Overwrite' : true|{false}. Overwrite previous tracking result.
%              'Frames' : Index array of frames to track (i.e., for downsampling). Default: all frames. 
%  'DownsamplingFactor' : Integer downsampling factor. Default: none.
%
%
% Example: runTracking(data, loadTrackSettings('Radius', [3 6], 'MaxGapLength', 2), 'Overwrite', true) ;

% Francois Aguet, May 2010 (last modified 05/28/2013)

function [] = runTracking(data, varargin)

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addOptional('settings', [], @isstruct);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('FileName', 'trackedFeatures', @ischar); % default of the tracker
ip.addParamValue('Frames', [], @isvector);
ip.addParamValue('DownsamplingFactor', [], @isscalar);
ip.addParamValue('DetectionFile', 'detection_v2.mat', @ischar);
ip.parse(data, varargin{:});
overwrite = ip.Results.Overwrite;
fileName = ip.Results.FileName;
idx = regexpi(fileName, '.mat');
if ~isempty(idx)
    fileName = fileName(1:idx-1);
end
frames = ip.Results.Frames;
dsfactor = ip.Results.DownsamplingFactor;
detectionFile = ip.Results.DetectionFile;

% Determine file name
if ~isempty(frames)
    fileName = [fileName '_customFrames(' num2str(frames(1)) '_' num2str(frames(end)) ')'];
end
if isempty(frames) && ~isempty(dsfactor)
    frames = 1:dsfactor:data.movieLength;
    fileName = [fileName '_' num2str(data.framerate*dsfactor) 's'];
end
fileName = [fileName '.mat'];

% Load tracker settings
settings = ip.Results.settings;
if isempty(settings)
   settings = loadTrackSettings('Radius', [3 6], 'MaxGapLength', 2);
end

% Run tracker on each data set
parfor i = 1:length(data)
    if ~(exist([data(i).source 'Tracking'], 'dir')==7) || overwrite
        fprintf('Running tracker on %s\n', getShortPath(data(i)));
        main(data(i), settings, fileName, detectionFile, frames);
    else
        fprintf('Tracking has already been run for %s\n', getShortPath(data(i)));
    end
end


function main(data, settings, fileName, detectionFile, frames)

dfile = [data.source 'Detection' filesep detectionFile];
if exist(dfile, 'file')==2
    dfile = load(dfile);
    movieInfo = dfile.frameInfo;
    if ~isempty(frames)
        movieInfo = movieInfo(frames);
    end
else
    fprintf(['runTracking: no detection data found for ' getShortPath(data)]);
    return;
end

[~,~] = mkdir([data.source 'Tracking']);
saveResults.dir = [data.source 'Tracking' filesep];
if ~isempty(fileName)
    saveResults.filename = fileName;
end
trackCloseGapsKalmanSparse(movieInfo, settings.costMatrices, settings.gapCloseParam, settings.kalmanFunctions, 2, saveResults, 1);
