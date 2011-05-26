function [] = runTracking(data, settings, varargin)
% runTracking tracks movies under a given condition folder. This creates
% the TrackInfo and lftInfo data structures necessary for lifetime
% analysis.
%
% SYNOPSIS [data] = runTracking(data,tracksettings,overwrite)
%
% INPUT     data(optional): structure containing the field source,
%               which specifies the directory for each movie (default is to
%               ask you to load them via gui)
%           tracksettings(optional):
%           overwrite(optional): 1 to overwrite to rerun tracking, 0 otherwise
%               (default is 0)
%
%          {'Frames', f} : array of frame indexes to track
%          {'DownsamplingFactor', d} : downsampling factor, must be an integer
%
% OUTPUT
%
% REMARKS   The function only performs the tracking for a given movie if
%           the folder Detection is not already present in the
%           specified directory; if you want a partial or faulty tracking
%           to be replaced, you need to DELETE these folders first.
%
%
% Francois Aguet, May 2010 (last modified 05/13/2011)

% Parse inputs
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
ip.addOptional('settings', []);
ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('FileName', 'trackedFeatures', @ischar); % default of the tracker
ip.addParamValue('Frames', [], @isvector);
ip.addParamValue('DownsamplingFactor', [], @isscalar);
ip.parse(data, varargin{:});
overwrite = ip.Results.Overwrite;
fileName = ip.Results.FileName;
frames = ip.Results.Frames;
dsfactor = ip.Results.DownsamplingFactor;

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
if isempty(settings)
    %load track settings required for tracking
    [fileName filePath] = uigetfile('.mat', 'Choose track settings mat file for tracking');
    load([filePath fileName]);
elseif ischar(settings)
    load(settings);
end

% Run tracker on each data set
parfor i = 1:length(data)
    if ~(exist([data(i).source 'Tracking'], 'dir')==7) || overwrite
        fprintf('Running tracker on %s\n', getShortPath(data(i)));
        main(data(i), settings, fileName, frames);
    else
        fprintf('Tracking has already been run for %s\n', getShortPath(data(i)));
    end
end



function [data] = main(data, settings, fileName, frames)

% now we're missing the variable movieInfo, which is the detection
% data. If a valid detection structure is a field in data, read it,
% else load the structure from the appropriate detection file
if isfield(data, 'detection') && ~isempty(data.detection)
    movieInfo = data.detection;
else
    %loadfile = load([data.source 'Detection' filesep 'detectionResults.mat']);
    loadfile = load([data.source 'Detection' filesep 'detection_v2.mat']);
    if isfield(loadfile, 'frameInfo')
        movieInfo = loadfile.frameInfo;
        if ~isempty(frames)
           movieInfo = movieInfo(frames); 
        end
    else
        error('No detection data file of specified format found');
    end
end

[~,~] = mkdir([data.source 'Tracking']);
saveResults.dir = [data.source 'Tracking' filesep];
if ~isempty(fileName)
    saveResults.filename = fileName;
end
trackCloseGapsKalmanSparse(movieInfo, settings.costMatrices, settings.gapCloseParam, settings.kalmanFunctions, 2, saveResults, 1);
