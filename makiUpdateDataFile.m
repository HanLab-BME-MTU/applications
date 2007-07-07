function makiUpdateDataFile(update)
%MAKIUPDATEDATAFILE is a developer tool to update maki data structures
%
% SYNOPSIS: makiUpdateDataFile(update)
%
% INPUT update: type of update
%           1: add cropping
%           2: reset cropping
%           3: remove faulty underscore in initCoordName
%           4: correct timeLapse in dataProperties
%           5: add new fields to dataStruct
%
% OUTPUT
%
% REMARKS This is YAH (yet another hack) that is appended whenever a new
%           update has to be run
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: jdorn
% DATE: 30-Jun-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test input
if nargin == 0 || isempty(update)
    disp('try again')
    return
end

% ask for top directory
topDir = uigetdir('','please select top directory');

% lookfor maki files
fileList = searchFiles('-makiData-','log',topDir,1);

% loop through fileList. Have update switch within loop to make life easier
nFiles = size(fileList,1);

for iFile = 1:nFiles
    % load the data file
    dataStruct = makiLoadDataFile(fullfile(fileList{iFile,2},fileList{iFile,1}));

    % remove all the xxx/xxxName pair files so that we don't get into
    % trouble with secureSave
    fn = fieldnames(dataStruct);
    for i=1:length(fn)
        % load data files that exist. Rest will be empty
        fileName = [fn{i},'Name'];
        if any(strmatch(fileName,fn)) && ~isempty(dataStruct.(fn{i}))
            try
                delete(fullfile(dataStruct.dataFilePath,dataStruct.(fileName)));

            catch
                warning('couldn''t delete %s',...
                    fullfile(dataStruct.dataFilePath,dataStruct.(fileName))); %#ok<WNTAG>
            end
        end
    end

    %===========================
    %----- UPDATE --------------
    %===========================

    % switch according to update
    switch update
        case 1
            % add cropImage as first entry in status
            dataStruct.status = [0; dataStruct.status];
            dataStruct.statusHelp = [{'cropMovie',''};dataStruct.statusHelp];

        case 2
            % reset cropping
            dataStruct.status(1) = -1;

        case 3
            dataStruct.initCoordName(end-4) = '';

        case 4
            % correct timeLapse in dataProperties
            dataStruct.dataProperties = ...
                defaultDataProperties(dataStruct.dataProperties);

        case 5

            % add fields
            if ~isfield(dataStruct,'tracks')
                projectName = dataStruct.projectName;
                dataStruct.tracksName=['tracks_',projectName,'.mat'];
                dataStruct.tracks=[];
                dataStruct.planeFitName=['planeFit_',projectName,'.mat'];
                dataStruct.planeFit=[];
                dataStruct.sisterListName=['sisterList_',projectName,'.mat'];
                dataStruct.sisterList=[];
                dataStruct.status = [dataStruct.status;[0,0,0]'];
                dataStruct.statusHelp = [dataStruct.statusHelp;cell(3,3)];
            end

        otherwise
            % do nothing
    end

    % save dataFile
    makiSaveDataFile(dataStruct);
end