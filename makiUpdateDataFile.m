function makiUpdateDataFile(update)
%MAKIUPDATEDATAFILE is a developer tool to update maki data structures
%
% SYNOPSIS: makiUpdateDataFile(update)
%
% INPUT update: type of update
%           1: add cropping.
%           2: reset cropping.
%           3: remove faulty underscore in initCoordName.
%           4: correct timeLapse in dataProperties.
%           5: add new fields to dataStruct.
%           6: swap track/planeFit in .statusHelp.
%           7: add defaults for groupSisters.
%           8: add _1 to data filenames in dataStruct.
%           9: update dataStruct to the actually lastest files in the
%              directory.
%           10: add .history field to dataStruct in the data file.
%           11: The old version of path 9 introduced wrong filenames into
%               the data structure if the file does not exist; this generated
%               files with wrong names in the directory. This patch fixes both,
%               the filenames in the dataStruct __and__ renames wrong files in
%               the working directory.
%           12: Swap entries in status, statusHelp and history to put slist
%               as 4th entry. Also add 8th entry for frameAlignment.
%           13: Add 2 new fields to dataStruct for frame alignment and
%               rearrange sequence of fields in dataStruct.
%           14: remove field goodTrackRatio from
%               dataStruct.dataProperties.groupSisters (basically I'm
%               fixing a bug).
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

if isempty(topDir)
    warning('No valid top directory selected');
    return;
end

% lookfor maki files
fileList = searchFiles('-makiData-','log',topDir,1);

% loop through fileList. Have update switch within loop to make life easier
nFiles = size(fileList,1);

for iFile = 1:nFiles
    % load the data file
    dataStruct = makiLoadDataFile(fullfile(fileList{iFile,2},fileList{iFile,1}));


    if update >= 9
        % remove all the xxx/xxxName pair files so that we don't get into
        % trouble with secureSave
        %
        % update 9 relies on the original files to adjust the broken data
        % structure -- thus do not delete
        %
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

        case 6
            % swap statusHelp-entries
            if strcmp(dataStruct.statusHelp{5,1},'tracks')
                dataStruct.statusHelp = ...
                    dataStruct.statusHelp([1,2,3,4,6,5,7,8],:);
            end

        case 7
            % add groupSisters - properties
            groupSisters.maxDist = 4; % maximum average sister separation in um
            groupSisters.goodTrackRatio = 0.75; % minimum relative track length for grouping
            groupSisters.robust = false; % whether or not to use robust statistics for determining cost function parameters
            groupSisters.costFunction = 'metaphase'; % cost function type
            dataStruct.dataProperties.groupSisters = groupSisters;

        case 8
            dataFn = {'dataPropertiesName', 'initCoordName', 'tracksName', ...
                'planeFitName', 'sisterListName', 'slistName' };
            for i=1:length(dataFn)
                % get data filenames
                projectName = dataStruct.projectName;
                fieldName = dataStruct.(dataFn{i});
                beginProjectName = regexp(fieldName,projectName);
                beginExt = regexp(fieldName,'.mat');
                if (beginProjectName + length(projectName)) == beginExt
                    % the project name is the last token in the field name
                    % Therefore '_1' needs to be added
                    dataStruct.(dataFn{i}) = [fieldName(1:(beginProjectName + length(projectName)-1)),...
                        '_1','.mat'];
                end
            end


        case 9
            dataFn = {'dataPropertiesName', 'initCoordName', 'tracksName', ...
                'planeFitName', 'sisterListName', 'slistName' };
            olddir = cd(dataStruct.dataFilePath);
            for i = 1:length(dataFn)
                searchString = dataFn{i};
                searchString = searchString(1:end-4);
                fileList2 = searchFiles(searchString);
                if ~isempty(fileList2)
                    % there are data files in the directory
                    fieldVersions = [];
                    for j = 1:size(fileList2,1)
                        fieldVersionTmp = makiGetVersion(fileList2{j,1});
                        if ~isempty(fieldVersionTmp)
                            % this is a file with a correct version tag
                            fieldVersions = [fieldVersions, fieldVersionTmp];
                        else
                            % this is a file with no version tag
                            fieldVersions = [fieldVersions, 0];
                        end
                    end
                    [maxVers,maxVersIndx]=max(fieldVersions);
                    if maxVers > 0
                        % adjust field names and load the corresponding
                        % file into the data structure
                        dataStruct.(dataFn{i}) = fileList2{maxVersIndx,1};
                        % fName = dataStruct.(dataFn{i});
                        % tmp = load(fullfile(dataStruct.dataFilePath,fName{1}));
                        tmp = load(fullfile(dataStruct.dataFilePath,dataStruct.(dataFn{i})));
                        fnTmp = fieldnames(tmp);
                        dataStruct.(searchString) = tmp.(fnTmp{1});
                    else
                        error(sprintf('%s is missing a version token. Use makiUpdateDataFile(8) first', searchString));
                    end
                else
                    % there is no such data in the directory yet
                    dataStruct.(dataFn{i}) = [searchString, '_', dataStruct.projectName, '_1', '.mat'];
                    dataStruct.(searchString) = [];
                end
            end
            cd(olddir);

        case 10
            if ~isfield(dataStruct,'history')
                dataFn = {'dataPropertiesName', 'initCoordName', 'planeFitName', ...
                    'tracksName', 'sisterListName', 'slistName' };

                if ~any(dataStruct.status(3:end)==1)
                    % nothing has been run yet on this movie (except
                    % cropping)
                    dataStruct.history = struct('numRuns',0,...
                        'dataProperties',[],...
                        'initCoord',[],...
                        'planeFit',[],...
                        'tracks',[],...
                        'sisterList',[],...
                        'slist',[]);
                else
                    % get newest version for each data field
                    dataStruct.history.numRuns = 1;
                    for i=1:length(dataFn)
                        % is there a data file for this field?
                        fieldDataName = dataFn{i};
                        fieldDataName = fieldDataName(1:end-4);
                        if ~isempty(dataStruct.(fieldDataName))
                            % get data filenames
                            fieldFileName = dataStruct.(dataFn{i});
                            fieldVersion = makiGetVersion(fieldFileName);
                            if ~isempty(fieldVersion)
                                dataStruct.history.(fieldDataName)=fieldVersion;
                            else
                                error(sprintf('%s is missing a version token. Use makiUpdateDataFile(8) first', fieldName));
                            end
                        else
                            % this task has no result yet
                            dataStruct.history.(fieldDataName)=NaN;
                        end
                    end
                end
            end
            
        case 11
            dataFn = {'dataPropertiesName', 'initCoordName', 'tracksName', ...
                'planeFitName', 'sisterListName', 'slistName' };
            for i=1:length(dataFn)
                if strmatch(dataFn{i},dataStruct.(dataFn{i}))
                    % the data filename wrongly contains the data field
                    % name
                    fileList2 = searchFiles(dataFn{i},'',dataStruct.dataFilePath);
                    fileNameBody = dataFn{i};
                    fileNameBody = fileNameBody(1:end-4);
                    if ~isempty(fileList2)
                        % one or more files in the data directory have the
                        % wrong filename
                        for j = 1:size(fileList2,1)
                            % rename files in directory
                            correctFileName = fileList2{j,1};
                            correctFileName = [fileNameBody,'_',correctFileName(length(dataFn{i})+1:end)];
                            movefile(fullfile(dataStruct.dataFilePath,fileList2{j,1}), ...
                                fullfile(dataStruct.dataFilePath,correctFileName));
                        end;
                    end
                    % update filenanme in data filename field independent
                    % of whether there was a file with the wrong filename
                    % in the directory
                    correctFileName = dataStruct.(dataFn{i});
                    correctFileName = [fileNameBody,'_',correctFileName(length(dataFn{i})+1:end)];
                    dataStruct.(dataFn{i}) = correctFileName;
                end
            end
            
        case 12
            %swap entries in status, statusHelp and history (entry 7 
            %becomes 4, shift the rest by 1) and add 8th entry
            
            %status
            tmpVar = dataStruct.status;
            if length(tmpVar) == 7
                tmpVar = [tmpVar(1:3); tmpVar(7); tmpVar(4:6); 0];
                dataStruct.status = tmpVar;
            end
            clear tmpVar
            
            %statusHelp
            tmpVar = dataStruct.statusHelp;
            if size(tmpVar,1) == 7
                tmpVar2 = cell(1,3);
                tmpVar2{1,1} = 'frameAlignment';
                tmpVar = [tmpVar(1:3,:); tmpVar(7,:); tmpVar(4:6,:); tmpVar2];
                dataStruct.statusHelp = tmpVar;
                clear tmpVar2
            end
            clear tmpVar
            
            %history - swap and replace NaNs with 0
            tmpVar2 = dataStruct.history;
            if ~isempty(tmpVar2)
                tmpVar.numRuns = tmpVar2.numRuns; %numRuns
                tmpVar3 = tmpVar2.dataProperties; %dataProperties
                tmpVar3(isnan(tmpVar3)) = 0;
                tmpVar.dataProperties = tmpVar3;
                tmpVar3 = tmpVar2.initCoord; %initCoord
                tmpVar3(isnan(tmpVar3)) = 0;
                tmpVar.initCoord = tmpVar3;
                tmpVar3 = tmpVar2.slist; %slist
                tmpVar3(isnan(tmpVar3)) = 0;
                tmpVar.slist = tmpVar3;
                tmpVar3 = tmpVar2.planeFit; %planeFit
                tmpVar3(isnan(tmpVar3)) = 0;
                tmpVar.planeFit = tmpVar3;
                tmpVar3 = tmpVar2.tracks; %tracks
                tmpVar3(isnan(tmpVar3)) = 0;
                tmpVar.tracks = tmpVar3;
                tmpVar3 = tmpVar2.sisterList; %sisterList
                tmpVar3(isnan(tmpVar3)) = 0;
                tmpVar.sisterList = tmpVar3;
                tmpVar.frameAlignment = zeros(1,tmpVar.numRuns); %frameAlignment
                dataStruct.history = tmpVar;
                clear tmpVar tmpVar3
            end
            clear tmpVar2
            
        case 13
            %add 2 new fields to dataStruct for frame alignment
            projectName = dataStruct.projectName;
            if ~isfield(dataStruct,'frameAlignmentName')
                dataStruct.frameAlignmentName = ['frameAlignment_' projectName '_1.mat'];
                dataStruct.frameAlignment = [];
            end
            %rearrange sequence of fields in dataStruct
            tmpVar = dataStruct;
            clear dataStruct;
            dataStruct.projectName = tmpVar.projectName;
            dataStruct.rawMovieName = tmpVar.rawMovieName;
            dataStruct.rawMoviePath = tmpVar.rawMoviePath;
            dataStruct.dataFileName = tmpVar.dataFileName;
            dataStruct.dataFilePath = tmpVar.dataFilePath;
            dataStruct.dataPropertiesName = tmpVar.dataPropertiesName;
            dataStruct.dataProperties = tmpVar.dataProperties;
            dataStruct.movieHeaderName = tmpVar.movieHeaderName;
            dataStruct.movieHeader = tmpVar.movieHeader;
            dataStruct.initCoordName = tmpVar.initCoordName;
            dataStruct.initCoord = tmpVar.initCoord;
            dataStruct.slistName = tmpVar.slistName;
            dataStruct.slist = tmpVar.slist;
            dataStruct.planeFitName = tmpVar.planeFitName;
            dataStruct.planeFit = tmpVar.planeFit;
            dataStruct.tracksName = tmpVar.tracksName;
            dataStruct.tracks = tmpVar.tracks;
            dataStruct.sisterListName = tmpVar.sisterListName;
            dataStruct.sisterList = tmpVar.sisterList;
            dataStruct.frameAlignmentName = tmpVar.frameAlignmentName;
            dataStruct.frameAlignment = tmpVar.frameAlignment;
            dataStruct.status = tmpVar.status;
            dataStruct.statusHelp = tmpVar.statusHelp;
            dataStruct.history = tmpVar.history;
            clear tmpVar
            
        case 14
            %remove dataStruct.dataProperties.groupSisters.goodTrackRatio
            
            dataProperties = dataStruct.dataProperties;
            if isfield(dataProperties,'groupSisters')
                groupSisters = dataProperties.groupSisters;
                if isfield(groupSisters,'goodTrackRatio')
                    groupSisters = rmfield(groupSisters,'goodTrackRatio');
                    dataStruct.dataProperties.groupSisters = groupSisters;
                end
            end
                
        otherwise
            % do nothing
    end

    % save dataFile
    makiSaveDataFile(dataStruct);
end
            