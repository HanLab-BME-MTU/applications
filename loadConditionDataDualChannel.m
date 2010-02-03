function [experiment] = loadConditionDataDualChannel(condDir)
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
% Francois Aguet, 12/14/2009

% select directory where all data for this condition are located
if (nargin == 0)
    condDir = [uigetdir(pwd, 'Select the ''condition'' folder') filesep];
end;
fprintf('Condition selected: %s\n', condDir);
    
% get directories where all data for this condition are located
expDir = dirList(condDir);

ct = 1;
if ~isempty(expDir)
    
    % get all experiment/date directories for this condition
    nExp = length(expDir);
    dateList(1:nExp) = struct('date', [], 'name', []);
    for i = 1:nExp
        
        % extract date from file name
        currDate = regexp(expDir(i).name, '\d+', 'match');
        if isempty(currDate)
            currDate = {'010101'}; % arbitrary default
        elseif (length(currDate) > 1)
            lengths = cellfun(@length, currDate);
            currDate = currDate(lengths == max(lengths(:)));
        end;
        dateList(i).date = currDate{1};
        
        % look for the individual cell data in this folder
        expPath = [condDir expDir(i).name filesep];
        cellDir = dirList(expPath);
        
        % loop over all cell folders
        if ~isempty(cellDir)
            for k = 1:length(cellDir)
                
                % extract framerate from directory name. Default to slow                
                if ( ~isempty(findstr(cellDir(k).name, 'fast')) || (~isempty(findstr(cellDir(k).name, '400ms'))) )
                    framerate = 0.4;
                else
                    framerate = 2;
                end
                
                if (ct == 1)
                    % get directories for individual channels
                    cellPath = [expPath cellDir(k).name filesep];
                    channel1Path = [uigetdir(cellPath, 'select first (master) channel (e.g. CCP channel)') filesep];
                    channel2Path = [uigetdir(cellPath, 'select second (slave) channel (e.g. other protein)') filesep];
                    channel1Name = channel1Path(length(cellPath)+1:end-1);
                    channel2Name = channel2Path(length(cellPath)+1:end-1);
                    fprintf('Channel 1 name: "%s"\n', channel1Name);
                    fprintf('Channel 2 name: "%s"\n', channel2Name);
                else
                    
                    % look for the individual channels in cell folder
                    if ~isempty(channel1Name)
                        channel1Path = [expPath cellDir(k).name filesep channel1Name filesep];
                    else
                        channel1Path = [expPath cellDir(k).name filesep];
                    end
                    if ~(exist(channel1Path, 'dir')==7)
                        %channel1Path = [uigetdir(cellDir(k).name, 'select first (master) channel (e.g. CCP channel)') filesep];
                    end
                    channel2Path = [expPath cellDir(k).name filesep channel2Name filesep];
                    if ~(exist(channel2Path, 'dir')==7)
                        %channel2Path = [uigetdir(cellDir(k).name, 'select second (slave) channel (e.g. other protein)') filesep];
                    end
                end
                
                fprintf('Loading: %s\n', channel1Path);

                experiment(ct).source = channel1Path;
                experiment(ct).channel1 = channel1Path;
                experiment(ct).channel2 = channel2Path;
                experiment(ct).date = currDate{1};
                experiment(ct).framerate = framerate;
                tifFiles = dir([experiment(ct).channel1 '*.tif']);
                experiment(ct).imagesize = size(imread([experiment(ct).channel1 tifFiles(1).name]));
                experiment(ct).movieLength = length(tifFiles);

                ct = ct+1;
            end;
        end;
    end;
else
    error('No data found in directory.');
end;