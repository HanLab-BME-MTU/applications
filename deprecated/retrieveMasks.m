function [maskArray] = retrieveMasks(experiment)

% retrieveMasks gets cell masks from CellAreaMasks folders for each
% experiment
%
% Input:
%           experiment=   data structure pointing to all the
%                       image/detection/tracking data for a given condition;
%                       for e.g. 10 cell movies, the structure should have 10
%                       entries; the data structure can be created with the
%                       function loadConditionData, and it needs to contain
%                       at least the field
%                       .source, which is the (path) location of the
%                       lifetime information folder
%                       .framerate, which is the movie framerate, which is
%                       necessary for the lifetime restriction
%
% OUTPUT
%           maskArray:  contains masks for each movie as a matrix of 3
%           dimensions where the third dimension represents the time axis;
%           of more than one movie, the output is a cell array where each
%           cell represents a different movie
% Uses:
%
% Daniel Nunez, updated May 06, 2009

%for each movie
for iexp = 1:length(experiment)

    clear masks

    %LOAD MOVIE DATA
    %go to cell directory
    cd(experiment(iexp).source)
    %load movie data
    load movieData.mat

    %LOAD MASKS
    %cd to mask directory
    cd(movieData.masks.directory)
    %find all tiffs
    files = dir('*.tif');
    %load mask
    for ifile = 1:length(files)
        masks(:,:,ifile) = imread(files(ifile).name);
    end

    %store results in array if more than one movie
    if length(experiment) == 1
        maskArray{iexp} = masks;
        %store results in matrix if one movie
    else
        maskArray{iexp} = masks;
    end

end %of for each experiment

end %of function