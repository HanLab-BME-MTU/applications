function dataStruct = makiLoadDataFile(dataFileName)
%MAKILOADDATAFILE loads a maki-data file from disk
%
% SYNOPSIS: dataStruct = makiLoadDataFile(dataFileName)
%
% INPUT dataFileName (opt): [pathName,filesep,filename] of dataFile. If
%           omitted, or if only a pathName is specified, the file can be
%           selected interactively
%
% OUTPUT dataStruct: data structure as described in makiMakeDataStruct
%
% REMARKS
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: jdorn
% DATE: 28-Jun-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%=================
%% CHECK INPUT
%=================

guiLoad = true;
guiPath = pwd;

if nargin < 1 || isempty(dataFileName)
    % guiLoad
else
    % check for file, then for path, then exit with error 
    if exist(dataFileName,'file')
        % we have a dataFile or a path
        [dataFilePath,part1,part2] = fileparts(dataFileName);
        if ~isempty(part2)
            % it's a data file (allow for stuff like .mat_old)
        dataFileName = [part1,part2];
        guiLoad = false;
        else
            % it's a path. fileParts will consider the last directory as
            % fileName
            guiPath = fullfile(dataFilePath,part1);
        end
    else 
        error('file or path %s not found',dataFileName);
    end
    
end


%=================
%% GUILOAD
%=================

if guiLoad
    % there is a bug in Matlab 2007a that doesn't allow the use of filter
    % names
%     oldDir = cd(guiPath);
%      [dataFileName, dataFilePath] = uigetfile(...
%          {'*-makiData-*','dataFiles';'*','all files'},'Please select makiDataFile');
%      cd(oldDir);
oldDir = cd(guiPath);
     [dataFileName, dataFilePath] = uigetfile(...
         '*-makiData-*','Please select makiDataFile');
     cd(oldDir);

    if dataFileName == 0
        error('no dataFile loaded')
    end
end

%================
%% LOAD DATA
%================

% load data file
load(fullfile(dataFilePath,dataFileName));
% write new path
dataStruct.dataFilePath = dataFilePath;
dataStruct.dataFileName = dataFileName;

% interpret rawMoviePath
if isempty(dataStruct.rawMoviePath)
    dataStruct.rawMoviePath = dataStruct.dataFilePath;
else
    % remove identifier
    dataStruct.rawMoviePath = makiPathDef(dataStruct.rawMoviePath);
end
        

% loop through dataStruct and read individual files
fn = fieldnames(dataStruct);
for i=1:length(fn)
    % load data files that exist. Rest will be empty
    fileName = [fn{i},'Name'];
    if any(strmatch(fileName,fn))
        try
            tmp = load(fullfile(dataStruct.dataFilePath,dataStruct.(fileName)));
            fnTmp = fieldnames(tmp);
            dataStruct.(fn{i}) = tmp.(fnTmp{1});
        catch
            dataStruct.(fn{i}) = [];
        end
    end
end
