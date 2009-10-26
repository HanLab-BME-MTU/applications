function [] = runTracking(experiment,tracksettings,force)
% runTracking tracks movies under a given condition folder. This creates
% the TrackInfo and lftInfo data structures necessary for lifetime
% analysis.
%
% SYNOPSIS [experiment] = runTracking()
%
% INPUT
%
% OUTPUT
%
% REMARKS   The function only performs the tracking for a given movie if
%           there is no tracking data (two folders called LifetimeInfo
%           and TrackInfoMatrices) already present in the specified
%           directory; if you want a partial or faulty tracking to be
%           replaced, you need to DELETE these folders first.
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

% CONVERT
% convert data from Henry's format to a format that is taken by the tracker
% if this is not necessary (a detection.mat file exists under a folder
% DetectionStructures) the function loadAndSaveDetection will not do
% anything to that movie
[experiment] = loadAndSaveDetection(experiment);

% add track settings to each experiment
for iexp = 1:length(experiment)
    experiment(iexp).tracksettings = tracksettings;
end

% TRACK
% this function creates the TrackInfo and lftInfo data structures if
% necessary
trackMissingFields(experiment,force);

end %of function