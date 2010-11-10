function [] = newTracker(data, settings, overwrite, fileName, frames)
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
% OUTPUT
%
% REMARKS   The function only performs the tracking for a given movie if
%           the folder Detection is not already present in the
%           specified directory; if you want a partial or faulty tracking
%           to be replaced, you need to DELETE these folders first.
%
% Dependencies  loadIndividualMovies
%               trackMissingFields
%
% Francois Aguet, May 2010


if nargin < 1 || isempty(data)
    % trim data structure as specified by user input
    [data] = loadIndividualMovies();
end
if nargin < 2 || isempty(settings)
    %load track settings required for tracking
    [fileName filePath] = uigetfile('.mat','choose track settings mat file for tracking');
    load([filePath fileName]);
elseif ischar(settings)
    load(settings);
end
if nargin < 3 || isempty(overwrite)
    overwrite = 0;
end
if nargin<4
    fileName = [];
end
if nargin<5
    frames = [];
end


parfor k = 1:length(data)
    fprintf('Tracking movie no. %d\n', k);
    if ~(exist([data(k).source 'Tracking'], 'dir')==7) || overwrite
        trackMissingFieldsNewTracker(data(k), settings, fileName, frames);
    else
        fprintf('Movie no. %d was skipped because it has already been tracked\n', k);
    end
end



function [data] = trackMissingFieldsNewTracker(data, settings, fileName, frames)

% now we're missing the variable movieInfo, which is the detection
% data. If a valid detection structure is a field in data, read it,
% else load the structure from the appropriate detection file
if isfield(data, 'detection') && ~isempty(data.detection)
    movieInfo = data.detection;
else
    loadfile = load([data.source 'Detection' filesep 'detectionResults.mat']);
    if isfield(loadfile, 'frameInfo')
        movieInfo = loadfile.frameInfo;
        if ~isempty(frames)
           movieInfo = movieInfo(frames); 
        end
    else
        error('No detection data file of specified format found');
    end
end
if ~(exist([data.source 'Tracking'], 'dir')==7)
    mkdir([data.source 'Tracking']);
end;
saveResults.dir = [data.source 'Tracking' filesep];
if ~isempty(fileName)
    saveResults.filename = fileName;
end
trackCloseGapsKalmanSparse(movieInfo, settings.costMatrices, settings.gapCloseParam, settings.kalmanFunctions, 2, saveResults, 1);