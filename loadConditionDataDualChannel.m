function [experiment] = loadConditionDataDualChannel()
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

% get directories where all data for this condition are located
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
        
        % loop over all cell folders
        if ~isempty(cellDir)
            for k = 1:length(cellDir)
                
                % extract framerate from the name of the cell folder
                % NOTE: if there's no specific identification for fast, default to slow
                
                if ( ~isempty(findstr(cellDir(k).name, 'fast')) || (~isempty(findstr(cellDir(k).name, '400ms'))) )
                    currFramerate = 0.4;
                else
                    currFramerate = 2;
                end
                
                if (k == 1)
                    % get directories for individual channels
                    channel1Path = uigetdir(cellDir(k).name, 'select first (master) channel (e.g. CCP channel)');
                    channel2Path = uigetdir(cellDir(k).name, 'select second (slave) channel (e.g. other protein)');
                    channel1Dir = channel1Path(find(channel1Path==filesep, 1, 'last')+1:end);
                    channel2Dir = channel2Path(find(channel2Path==filesep, 1, 'last')+1:end);
                else
                    
                    % look for the individual channels in cell folder
                    if exist(channel1Dir, 'dir')==7
                        channel1Path = [cellDir(k).name filesep channel1Dir];
                    else
                        channel1Path = uigetdir(cellDir(k).name, 'select first (master) channel (e.g. CCP channel)');
                    end
                    
                    if exist(channel2Dir, 'dir')==7
                        channel2Path = [cellDir(k).name filesep channel2Dir];
                    else
                        channel2Path = uigetdir(cellDir(k).name, 'select second (slave) channel (e.g. other protein)');
                    end
                end
                
                % enter data
                experiment(ct).source =   channel1Path;
                experiment(ct).channel1 = channel1Path;
                experiment(ct).channel2 = channel2Path;
                experiment(ct).date = currDate;
                experiment(ct).framerate = currFramerate;
                
                ct = ct+1;
            end
        end
    end
else
    error('no usable data in directory');
end