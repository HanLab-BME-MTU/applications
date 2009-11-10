function [] = runTracking(experiment,tracksettings,force)
% runTracking tracks movies under a given condition folder. This creates
% the TrackInfo and lftInfo data structures necessary for lifetime
% analysis.
%
% SYNOPSIS [experiment] = runTracking(experiment,tracksettings,force)
%
% INPUT     experiment(optional): structure containing the field source,
%               which specifies the directory for each movie (default is to
%               ask you to load them via gui)
%           tracksettings(optional):
%           force(optional): 1 to force to rerun tracking, 0 otherwise
%               (default is 0)
%
% OUTPUT
%
% REMARKS   The function only performs the tracking for a given movie if
%           the folder DetectionStructures is not already present in the
%           specified directory; if you want a partial or faulty tracking
%           to be replaced, you need to DELETE these folders first.
%
% Dependencies  loadIndividualMovies
%               loadAndSaveDetection
%               trackMissingFields
%
% Daniel Nunez, March 5, 2008

%INPUTS
if nargin < 3 || isempty(force)
    force = 0;
end
if nargin < 1 || isempty(experiment)
    % trim experiment structure as specified by user input
    [experiment] = loadIndividualMovies();
end
if nargin < 2 || isempty(tracksettings)
    %load track settings required for tracking
    [fileName filePath] = uigetfile('.mat','choose track settings mat file for tracking');
    load([filePath filesep fileName]);
end

%ONLY TRACK IF DETECTIONSTRUCTURES DOES NOT EXIST. OTHERWISE MUST FORCE TO
%GET TO TRACK
parfor iexp = 1:length(experiment)
    %select experiment of interest
    exp = experiment(iexp);
    %go to experiment folder
    cd(exp.source);
    if ~ ( exist('DetectionStructures')==7 ) || force
        % CONVERT
        % convert data from Henry's format to a format that is taken by the tracker
        % if this is not necessary (a detection.mat file exists under a folder
        % DetectionStructures) the function loadAndSaveDetection will not do
        % anything to that movie
        [exp] = loadAndSaveDetection(exp);
        
        % add track settings to each experiment
        exp.tracksettings = tracksettings;
        
        % TRACK
        % this function creates the TrackInfo and lftInfo data structures if
        % necessary
        trackMissingFields(exp,force);
        
    end %of if DetectionSturctures exists
end %of for each experiment
end %of function