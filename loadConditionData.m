function [experiment] = loadConditionData(condDir, marker)
% loadConditionData loads the relevant information for all the data
% available for a specific experiment condition; this requires a specific
% dircetory structure and nomenclature (see below)
%
% SYNOPSIS [experiment] = loadConditionData()
%
% INPUT
%
% OUTPUT   experiment: structure with the fields
%                       .source = pathname of the data/movie
%                       .date = date when movie was taken
%                       .framerate = framerate of the movie (2s or 0.4s)
%
%
%
% Dinah Loerke, January 24th, 2008
% Francois Aguet, 11/02/2009


if nargin<1
    % select directory where all data for this condition are located
    condDir = [uigetdir(pwd, 'Please select the folder for this condition') filesep];
end

% get directories for experiments under this condition
expDir = dirList(condDir);

experiment = struct('source', [], 'channel1', [], 'date', [], 'framerate', [],...
    'imagesize', [], 'movieLength', [], 'channel1marker', []);

ct = 1;
if ~isempty(expDir)
    for i = 1:length(expDir)
        
        % extract date from file name
        currDate = regexp(expDir(i).name, '\d+', 'match');
        if isempty(currDate)
            currDate = {'010101'}; % arbitrary default
        elseif (length(currDate) > 1)
            lengths = cellfun(@length, currDate);
            currDate = currDate(lengths == max(lengths(:)));
        end;
        currDate = currDate{1};
        expPath = [condDir expDir(i).name];
        
        % look for the individual cell data in this folder
        cellDir = dirList(expPath);
        
        % loop over all cells
        if ~isempty(cellDir)
            for k = 1:length(cellDir)
                
                % extract framerate from the name of the cell folder
                % NOTE: if there's no specific identification for fast, default to slow
                if ~isempty(regexp(cellDir(k).name, '\d+s', 'match'))
                    fr = regexp(cellDir(k).name, '\d+s', 'match');
                    framerate = str2double(fr{1}(1:end-1));
                elseif ~isempty(regexp(cellDir(k).name, '\d+ms', 'match'))
                    fr = regexp(cellDir(k).name, '\d+ms', 'match');
                    framerate = str2double(fr{1}(1:end-2))/1000;
                elseif ~isempty(findstr(cellDir(k).name, 'fast'))
                    framerate = 0.4;
                else
                    framerate = 2;
                end
                
                % enter data
                experiment(ct).source = [expPath filesep cellDir(k).name filesep];
                experiment(ct).channel1 = experiment(ct).source;
                experiment(ct).date = currDate;
                experiment(ct).framerate = framerate;
                if nargin>1
                    experiment(ct).channel1marker = marker;
                end
                if nargin<3
                    experiment(ct).NA = 1.49;
                    experiment(ct).M = 100;
                    experiment(ct).pixelSize = 6.7e-6;
                end
                
                tifFiles = dir([experiment(ct).source '*.tif*']);
                if ~isempty(tifFiles)
                    experiment(ct).imagesize = size(imread([experiment(ct).source tifFiles(1).name]));
                    experiment(ct).movieLength = length(tifFiles);
                    ct = ct+1;
                end
            end
        end
    end
else
    error('no usable data in directory');
end