function updateProjects(varargin)
%UPDATEEXPERIMENTALDATA is a tool to change parts of or information in chromdyn experimental data structures
%
% SYNOPSIS updateExperimentalData(opt1,opt2,...)
%
% INPUT    opt1,opt2,... : strings telling what to update
%                           - 'changePath' : updates relative path of
%                                  project. This will be done, too, if any
%                                  of the two options below are selected
%                           - 'changePixelsize' : update from old to new
%                                  pixelsize (hardcoded change)
%                           - 'cleanupFiles' : retains only the most recent
%                                  slists/idlists
%
% The program will ask the user to select the directory the code should run
% on. This is a further check to make sure the person knows what's going on
%
% Please run this code with care - you could screw up all your work with
% one click
%
%c: jonas, 04/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%========================
% SET DEFAULTS
%========================

changePath = 0;
changePixelsize = 0;
pixelsizeOld = 0.000;
pixelsizeNew = 0.000;
cleanupFiles = 0;

%========================
% TEST INPUT
%========================

% check first that we have at least one input argument, then change the
% defaults. Remember what the user chose for the dialog, and exit if
% nothing is to be done

if nargin == 0
    error('please specify at least one input argument for UPDATEPROJECTS!')
end

% init actionString
actionString = '';

if any(strmatch('changePath',varargin,'exact'))
    changePath = 1;
    actionString = [actionString,' & update paths'];
end

if any(strmatch('changePixelsize',varargin,'exact'))
    changePixelsize = 1;
    actionString = [actionString,' & change pixelsizes'];
end

if any(strmatch('cleanupFiles',varargin,'exact'))
    cleanupFiles = 1;
    actionString = [actionString,' & cleanup old files'];
end

if isempty(actionString)
    % nothing matched
    error('no action selected in UPDATEPROJECTS. Please check your spelling!')
else
    % remove the '&'
    actionString = actionString(3:end);
    
    % add changePath as freebie (you can't break anything there)
    if (changePixelsize | cleanupFiles) & ~changePath
        changePath = 1;
        actionString = [actionString ' (includes updating path)'];
    end
end

%=====================
% ASK FOR INPUT
%=====================

selection = myQuestdlg(...
    ['You are about to',...
        actionString,...
        '. If this is not what you want to do, press ''Cancel''. ',...
        'Else, please select which directory (including subdirectories) will be changed.'],...
    'Please read carefully',...
    'Select Top Directory',...
    'Cancel',...
    'Cancel');

% test selection and quit if cancel
if strcmp(selection,'Select Top Directory')
    topDirectoryName = uigetdir;
    if topDirectoryName == 0
        error('directory selection cancelled by user')
    end
else
    error('directory selection cancelled by user')
end

%=======================================
% BUILD FILE/DIR-LIST, GET BIODATA-STUFF
%=======================================

listOfDataFiles = searchFiles('-data-','log',[],1,'new');
mainDir = cdBiodata(4);
mainDirLength = length(mainDir);



%====================
% DO UPDATE
%====================

% loop through the listOfDataFiles. First cleanup files, then change
% pixelsize, then do the path update

for iProject = 1:size(listOfDataFiles,1)
    
    % init var to remember changed vars, so that we only save the necessary
    % vars of the data file
    changedVariables = cell(0);
    
    % load the current project file
    currentDir = [listOfDataFiles{iProject,2}, filesep];
    projectDataFile = load([currentDir,listOfDataFiles{iProject,1}]);
    
    %==========
    % CLEANUP
    %==========
    
    % For cleanup, delete old idlists from project folder, then save
    % them out of the projectDataFile. (if changePixelSize, we do not 
    % save idlists until the pixelsize has been changed)
    % Also remove old slists, intlists, older project data files, 
    % older intFigures, older filtered movies and tmpDataProperties.
    
    if cleanupFiles
        
        %-----delete old slists
        slistList = dir([currentDir 'slist*']);
        % find the newest date
        slistDates = cat(1,slistList.date)
        [dummy,dateIdx] = sort(slistDates);
        % the newest date will be at the end of the list, so delete all
        % others
        for i = 1:length(dateIdx-1)
            delete([currentDir slistList(dateIdx(i)).name]);
        end
        
        %-----delete idlists - save later
        delete([currentDir 'idlist*']);
        
        %-----delete intLists
        delete([currentDir 'intList*']);
        
        %-----delete old intFigures
        intFigureList = dir([currentDir 'intFigure*']);
        % find the newest date
        intFigureDates = cat(1,intFigureList.date)
        [dummy,dateIdx] = sort(intFigureDates);
        % the newest date will be at the end of the list, so delete all
        % others
        for i = 1:length(dateIdx-1)
            delete([currentDir intFigureList(dateIdx(i)).name]);
        end
        
        %-----delete old data files
        dataFileList = dir([currentDir '*-data-*']);
        for i = 1:length(dataFileList)
            if ~strmatch(listOfDataFiles{iProject,1},dataFileList(i).name)
                delete([currentDir dataFileList(i).name]);
            end
        end
        
        %-----delete old filtered movies
        filteredMovieList = dir([currentDir '*.fim']);
        filteredMovieDates = cat(1,filteredMovieList.date)
        [dummy,dateIdx] = sort(filteredMovieDates);
        % the newest date will be at the end of the list, so delete all
        % others
        for i = 1:length(dateIdx-1)
            delete([currentDir filteredMovieList(dateIdx(i)).name]);
        end
        
        %------delete tmpDataProperties
        delete([currentDir 'tmpDataProperties.mat']);
        
        %-------------------------------------------
        % save idlists out of projectDataFile
        
        if ~(changePixelsize & projectDataFile.dataProperties.PIXELSIZE_XY == pixelsizeOld)
            if isfield(projectDataFile,'idlist')
                idlist = projectDataFile.idlist;
                idname = [currentDir 'idlist-' idlist(1).stats.created{1}];
                save(idname,'idlist');
                clear idlist
                
                % only if there was idlist, there will be idlist_L
                if isfield(projectDataFile,'idlist_L')
                    idlist_L = projectDataFile.idlist_L;
                    idname = [currentDir 'idlist-L-' idlist_L(1).stats.created{1}];
                    save(idname,'idlist_L');
                    clear idlist_L
                    
                    % only if there was idlist_L, there will be idlisttrack
                    if isfield(projectDataFile,'idlisttrack')
                        idlisttrack = projectDataFile.idlisttrack;
                        idname = [currentDir 'idlisttrack-' idlisttrack(1).stats.created{1}];
                        save(idname,'idlisttrack');
                        clear idlisttrack
                        
                        % only if there was idlisttrack, there will be idlisttrack_L
                        if isfield(projectDataFile,'idlisttrack_L')
                            idlisttrack_L = projectDataFile.idlisttrack_L;
                            idname = [currentDir 'idlisttrack-L-' idlisttrack_L(1).stats.created{1}];
                            save(idname,'idlisttrack_L');
                            clear idlisttrack_L
                            
                        end
                    end
                end
            end
        end
        
        
        
    end % if cleanupFiles
    
    %==========
    % PIXELSIZE
    %==========
    
    % check data properties if the pixelsizes are wrong. If yes, update
    % dataProperties and all the idlists. If no cleanup, load all the
    % idlists and update (otherwise, just write to file). Then, update the
    % r3dHeader and tmpDataProperties
    
    if changePixelsize & projectDataFile.dataProperties.PIXELSIZE_XY == pixelsizeOld;
        
        pixelsizeRatio = pixelsizeNew/pixelsizeOld;
        
        
        % change pixelsize in dataProperties
        projectDataFile.dataProperties.PIXELSIZE_XY = pixelsizeNew;
        changedVariables{end+1} = 'dataProperties';
        
        % change pixelsize in idlists on disk. Do this before saving
        % already changed stuff
        idlistList = dir([currentDir 'idlist*']);
        for iList = 1:length(idlistList)
            % load idlist
            id = load([currentDir idlistList(iList).name]);
            % find variable name
            idname = char(fieldnames(id));
            % assign as idlist
            eval(['idlist = id.' idname ';']);
            % change pixelsize
            for t = 1:length(idlist)
                if ~isempty(idlist(t).linklist)
                    idlist(t).linklist(:,9:10) = idlist(t).linklist(:,9:10)*pixelsizeRatio;
                end
            end
            % assign back
            eval([idname '= idlist;']);
            % save to disk - use name idname
            eval(['save([currentDir idlistList(iList).name],' idname ');']);
            % clear vars
            clear('idlist', idname)
        end
        
        %----------------
        % change pixelsize in idlists in projectDataFile. We have to loop through every
        % timepoint, unfortunately. If we were doing a cleanup, we will
        % save the idlist externally.
        
        if isfield(projectDataFile,'idlist')
            idlist = projectDataFile.idlist;
            for t = 1:length(idlist)
                if ~isempty(idlist(t).linklist)
                    idlist(t).linklist(:,9:10) = idlist(t).linklist(:,9:10)*pixelsizeRatio;
                end
            end
            idlist(1).stats.status{end+1} = [nowString 'changed XY pixelsize'];
            projectDataFile.idlist = idlist;
            changedVariables{end+1} = 'idlist';
            
            if cleanupFiles
                idname = [currentDir 'idlist-' idlist(1).stats.created{1}];
                save(idname,'idlist');
            end
            clear idlist
            
            % only if there was idlist, there will be idlist_L
            if isfield(projectDataFile,'idlist_L')
                idlist_L = projectDataFile.idlist_L;
                
                for t = 1:length(idlist_L)
                    if ~isempty(idlist_L(t).linklist)
                        idlist_L(t).linklist(:,9:10) = idlist_L(t).linklist(:,9:10)*pixelsizeRatio;
                    end
                end
                idlist_L(1).stats.status{end+1} = [nowString 'changed XY pixelsize'];
                projectDataFile.idlist_L = idlist_L;
                changedVariables{end+1} = 'idlist_L';
                
                if cleanupFiles
                    idname = [currentDir 'idlist-L-' idlist_L(1).stats.created{1}];
                    save(idname,'idlist_L');
                end
                clear idlist_L
                
                % only if there was idlist_L, there will be idlisttrack
                if isfield(projectDataFile,'idlisttrack')
                    idlisttrack = projectDataFile.idlisttrack;
                    
                    for t = 1:length(idlisttrack)
                        if ~isempty(idlisttrack(t).linklist)
                            idlisttrack(t).linklist(:,9:10) = idlisttrack(t).linklist(:,9:10)*pixelsizeRatio;
                        end
                    end
                    idlisttrack(1).stats.status{end+1} = [nowString 'changed XY pixelsize'];
                    projectDataFile.idlisttrack = idlisttrack;
                    changedVariables{end+1} = 'idlisttrack';
                    
                    if cleanupFiles
                        idname = [currentDir 'idlisttrack-' idlisttrack(1).stats.created{1}];
                        save(idname,'idlisttrack');
                    end
                    clear idlisttrack
                    
                    % only if there was idlisttrack, there will be idlisttrack_L
                    if isfield(projectDataFile,'idlisttrack_L')
                        idlisttrack_L = projectDataFile.idlisttrack_L;
                        
                        for t = 1:length(idlisttrack_L)
                            if ~isempty(idlisttrack_L(t).linklist)
                                idlisttrack_L(t).linklist(:,9:10) = idlisttrack_L(t).linklist(:,9:10)*pixelsizeRatio;
                            end
                        end
                        idlisttrack_L(1).stats.status{end+1} = [nowString 'changed XY pixelsize'];
                        projectDataFile.idlisttrack_L = idlisttrack_L;
                        changedVariables{end+1} = 'idlisttrack_L';
                        
                        if cleanupFiles
                            idname = [currentDir 'idlisttrack-L-' idlisttrack_L(1).stats.created{1}];
                            save(idname,'idlisttrack_L');
                        end
                        clear idlisttrack_L
                        
                    end
                end
            end
        end
        
        % change pixelsize in r3dMovieHeader
        try
            load([currentDir 'r3dMovieHeader'])
            r3dMovieHeader.pixelX = pixelsizeNew;
            r3dMovieHeader.pixelY = pixelsizeNew;
            save([currentDir 'r3dMovieHeader'], 'r3dMovieHeader');
            clear('r3dMovieHeader')
        catch
            if ~strcmp(lasterr ['Error using ==> load' char(10) 'Unable to read file ' currentDir 'r3dMovieHeader: file does not exist.'])
                % we care about all errors except the nonexistent header
                rethrow(lasterror)
            end
        end
        
        % change pixelsize in tmpDataProperties
        try
            load([currentDir 'tmpDataProperties'])
            tmpDataProperties.PIXELSIZE_XY = pixelsizeNew;
            save([currentDir 'tmpDataProperties'], 'tmpDataProperties');
            clear('tmpDataProperties')
        catch
            if ~strcmp(lasterr ['Error using ==> load' char(10) 'Unable to read file ' currentDir 'tmpDataProperties: file does not exist.'])
                % we care about all errors except the nonexistent header
                rethrow(lasterror)
            end
        end
        
    end % if changePixelsize
    
    %==========
    % PATH
    %==========
    
    % just set the path in projectProperties to the current filePath
    
    if changePath
        
        %update file location - check whether to use relative or absolute path
        if strmatch(lower(mainDir),lower(currentDir))
            if ~strmatch(lower(mainDir),lower(currentDir),'exact')
                projectDataFile.projProperties.dataPath = currentDir(mainDirLength+2:end);
            end
        else
            projectDataFile.projProperties.dataPath = currentDir;
        end
        changedVariables{end+1} = 'projProperties';
        
    end % if changePath
    
    %================
    % save data file
    %================
    
    for i = 1:length(changedVariables)
        
        % read variable from structure
        eval([changedVariables{i} '= projectDataFile.' changedVariables{i} ';']);
        
        % save
        eval(['save([currentDir,listOfDataFiles{iProject,1}],' changedVariables{i} ', ''-append'' );']);
        
    end
    
end % for iProject = 1:size(listOfDataFiles,1)
