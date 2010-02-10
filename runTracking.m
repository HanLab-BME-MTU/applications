function [] = runTracking(exp, tracksettings, overwrite)
% runTracking tracks movies under a given condition folder. This creates
% the TrackInfo and lftInfo data structures necessary for lifetime
% analysis.
%
% SYNOPSIS [exp] = runTracking(exp,tracksettings,overwrite)
%
% INPUT     exp(optional): structure containing the field source,
%               which specifies the directory for each movie (default is to
%               ask you to load them via gui)
%           tracksettings(optional):
%           overwrite(optional): 1 to overwrite to rerun tracking, 0 otherwise
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
%               trackMissingFields
%
% Daniel Nunez, March 5, 2008
% Francois Aguet, Jan 2010

%INPUTS
if nargin < 3 || isempty(overwrite)
    overwrite = 0;
end
if nargin < 1 || isempty(exp)
    % trim exp structure as specified by user input
    [exp] = loadIndividualMovies();
end
if nargin < 2 || isempty(tracksettings)
    %load track settings required for tracking
    [fileName filePath] = uigetfile('.mat','choose track settings mat file for tracking');
    load([filePath filesep fileName]);
end

detectFlag = zeros(1,length(exp));
for i = 1:length(exp)
    detectionFile = [exp(i).source 'DetectionStructures' filesep 'detection.mat'];
    if (~exist(detectionFile, 'file')==2) || overwrite
        if (~exist([exp(i).source 'DetectionStructures'], 'dir')==7)
            mkdir([exp(i).source 'DetectionStructures']);
        end;
        fprintf('Converting detection data for movie no. %d\n', i);
        % Convert data from Henry's format to the format read by the tracker
        [detection] = convertDetectDataForTracking([exp(i).source filesep 'maxdata283']);
        save detectionFile detection;
        detectFlag(i) = 1;
    end
end
useExp = exp(detectFlag==1);

% two separate loops are used because 'parfor' does not support the 'save' function.
parfor i = 1:length(useExp)
    iexp = useExp(i);
    iexp.tracksettings = tracksettings;
    % this function creates the TrackInfo and lftInfo data structures
    trackMissingFields(iexp, overwrite);
end