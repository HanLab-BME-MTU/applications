function [experiment] = loadConditionData()
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


% select directory where all data for this condition are located
condDir = uigetdir(pwd, 'Please select the folder for this condition');

% get directories for experiments under this condition
expDir = dirList(condDir);

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
        expPath = [condDir filesep expDir(i).name];

        % look for the individual cell data in this folder
        cellDir = dirList(expPath);
        
        % loop over all cells
        if ~isempty(cellDir)
            for k = 1:length(cellDir)

                % extract framerate from the name of the cell folder
                % NOTE: if there's no specific identification for fast, default to slow
                if ( ~isempty(findstr(cellDir(k).name, 'fast')) || ~isempty(findstr(cellDir(k).name, '400ms')) )
                    currFramerate = 0.4;
                else 
                    currFramerate = 2;
                end
                
                % enter data
                experiment(ct).source = [expPath filesep cellDir(k).name];
                experiment(ct).date = currDate;
                experiment(ct).framerate = currFramerate;
                
                tifFiles = dir([experiment(ct).source filesep '*.tif']);
                experiment(ct).imagesize = size(imread([experiment(ct).source filesep tifFiles(1).name]));
                experiment(ct).movieLength = length(tifFiles);

                ct = ct+1;
            end
        end   
    end
else
    error('no usable data in directory');
end      