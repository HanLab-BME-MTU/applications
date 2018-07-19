function [] = runOldTracking(data, tracksettings, overwrite)
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
% Daniel Nunez, March 5, 2008
% Francois Aguet, April 2010, updated November 2010

if nargin < 3 || isempty(overwrite)
    overwrite = 0;
end
if nargin < 1 || isempty(data)
    % trim data structure as specified by user input
    [data] = loadIndividualMovies();
end
if nargin < 2 || isempty(tracksettings)
    %load track settings required for tracking
    [fileName filePath] = uigetfile('.mat','choose track settings mat file for tracking');
    load([filePath fileName]);
elseif ischar(tracksettings)
    load(tracksettings);
end

nd = length(data);

parfor k = 1:length(data)
    idata = data(k);
    idata.tracksettings = tracksettings;
    fprintf('Tracking movie #%d/%d: %s\n', k, nd, idata.source);
    
    trackMissingFields(idata, overwrite);
end



function [data] = trackMissingFields(data,overwrite)

if nargin < 2 || isempty(overwrite)
    overwrite = 0;
end

for k = 1:length(data)
    
    if ~(exist([data(k).source 'TrackInfoMatrices'], 'dir')==7) || overwrite
        
        costMatrices    = data(k).tracksettings.costMat;
        gapCloseParam   = data(k).tracksettings.gapClosePar;
        iterParam       = data(k).tracksettings.iterPar;
        
        % now we're missing the variable movieInfo, which is the detection
        % data. If a valid detection structure is a field in data, read it,
        % else load the structure from the appropriate detection file
        if isfield(data(k), 'detection') && ~isempty(data(k).detection)
            movieInfo = data(k).detection;
        else
            loadfile = load([data(k).source 'Detection' filesep 'detectionResults.mat']);
            if isfield(loadfile, 'frameInfo')
                movieInfo = loadfile.frameInfo;
            else
                error('No detection data file of specified format found');
            end
        end
        if ~(exist([data(k).source 'TrackInfoMatrices'], 'dir')==7)
            mkdir([data(k).source 'TrackInfoMatrices']);
        end;
        saveResults.dir = [data.source 'TrackInfoMatrices' filesep] ;
        trackWithGapClosing(movieInfo, costMatrices, 'getTrackStats', gapCloseParam, iterParam, saveResults);
    else
        fprintf('Movie no. %d was skipped because it has already been tracked\n', k);
    end
end
