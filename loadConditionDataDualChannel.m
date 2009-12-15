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
    condDir = uigetdir(pwd, 'Select the ''condition'' folder');
end;
fprintf('Condition selected: %s\n', condDir);
    
% get directories where all data for this condition are located
expDir = dirList(condDir);

ct = 1;
if ~isempty(expDir)
    
    % get all experiment/date directories for this condition
    nExp = length(expDir);
    dateList(1:nExp) = struct('date', [], 'name', []);
    %nCells = zeros(1,nExp);
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
%         dateList(i).name = expDir(i).name;
        
        % look for the individual cell data in this folder
        expPath = [condDir filesep expDir(i).name];
        cellDir = dirList(expPath);
        %nCells(i) = length(cellDir);
        
        % loop over all cell folders
        if ~isempty(cellDir)
            for k = 1:length(cellDir)
                cellPath = [expPath filesep cellDir(k).name];
                
                % extract framerate from directory name. Default to slow                
                if ( ~isempty(findstr(cellDir(k).name, 'fast')) || (~isempty(findstr(cellDir(k).name, '400ms'))) )
                    framerate = 0.4;
                else
                    framerate = 2;
                end
                
                if (ct == 1)
                    % get directories for individual channels
                    channel1Path = uigetdir(cellPath, 'select first (master) channel (e.g. CCP channel)');
                    channel2Path = uigetdir(cellPath, 'select second (slave) channel (e.g. other protein)');
                    channel1Name = channel1Path(find(channel1Path==filesep, 1, 'last')+1:end);
                    channel2Name = channel2Path(find(channel2Path==filesep, 1, 'last')+1:end);
                    fprintf('Channel 1 name: %s\n', channel1Name);
                    fprintf('Channel 2 name: %s\n', channel2Name);
                else
                    
                    % look for the individual channels in cell folder
                    channel1Path = [cellPath filesep channel1Name];
                    if ~exist(channel1Path, 'dir')==7
                        channel1Path = uigetdir(cellDir(k).name, 'select first (master) channel (e.g. CCP channel)');
                    end
                    channel2Path = [cellPath filesep channel2Name];
                    if ~exist(channel2Path, 'dir')==7
                        channel2Path = uigetdir(cellDir(k).name, 'select second (slave) channel (e.g. other protein)');
                    end
                end
                
                % enter data
                experiment(ct).source = channel1Path;
                experiment(ct).channel1 = channel1Path;
                experiment(ct).channel2 = channel2Path;
                experiment(ct).date = currDate{1};
                experiment(ct).framerate = framerate;
                
                tifFiles = dir([experiment(ct).channel1 filesep '*.tif']);
                experiment(ct).imagesize = size(imread([experiment(ct).channel1 filesep tifFiles(1).name]));
                experiment(ct).movieLength = length(tifFiles);

                ct = ct+1;
            end;
        end;
    end;
else
    error('No data found in directory.');
end;