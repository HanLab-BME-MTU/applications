function [idlist,dataProperties,projectProperties,slist,filteredMovie] = loadProjectData(fileName,pathName,idname,GUI,forceNew,noMovie)
%LOADPROJECTDATA loads experimental Chromdyn data from file
%
% SYNOPSIS [idlist,dataProperties,projectProperties,slist,filteredMovie] = loadProjectData(fileName,pathName,idname)
%
% INPUT    fileName  (opt) name of data file. if 1, the program will load
%                       GUI mode. If empty, the program will attempt to
%                       load a dataFile from the current or specified path
%          pathName  (opt) name of path for data file. If fileName, but no
%                       pathName is specified, currentDir is used
%          idname    (opt) name of idlist to be loaded. Can be "last"
%          GUI       (opt) whether to show different idlist choices in GUI
%          forceNew  (opt) if 1, loadProjectData will only look for idX2.
%                       If there aren't any, idlist will be -1.
%          noMovie   (opt) if 1, loadStruct for filtered movie will be
%                       returned. If 2, loadStruct for corr/raw will be
%                       returned. If 3, filteredMovie will be emtpy
%
% OUTPUT   idlist           user selected idlist if several possible
%                           idlist(1).stats.idname contains current idname
%          dataProperties
%          projectProperties
%          slist
%          filteredMovie
%
%c: jonas 05/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%===================
% TEST INPUT
%===================

if nargin < 4 || isempty(GUI)
    GUI = 0;
end

if nargin == 0
    loadData = 1;
elseif isempty(fileName)
    loadData = 2;
else
    if fileName == 1
        loadData = 1;
        GUI = 1;
    else
        loadData = 0;
    end
end

if nargin < 2 || isempty(pathName)
    % if fileName already contains the path, we don't need to add pathName
    % in the beginning
    if ~isempty(strfind(fileName,filesep))
        pathName = [];
    else
        pathName = pwd;
    end
end
% if pathName is empty, we don't need to add filesep, of course
if ~isempty(pathName) && ~strcmp(pathName(end),filesep)
    pathName = [pathName,filesep];
end

if nargin < 3 || isempty(idname)
    idname = [];
end

%input4: GUI

if nargin < 5 || isempty(forceNew)
    forceNew = 0;
end

if nargin < 6 || isempty(noMovie)
    noMovie = 0;
end


%=================


%========================
% FIND FILE IF NECESSARY
%========================

switch loadData
    case 0 % check that the file exists
        if ~exist(fileName,'file') && ~ exist([pathName,fileName],'file')
            if GUI
                h = errordlg('file not found');
                uiwait(h)
                return
            else
                error('file not found')
            end
        end

    case 1 % launch GUI
        oldDir = pwd;

        %check if default biodata-dir exists and cd if exist
        mainDir = cdBiodata(2);

        %get project data file
        [fileName,pathName] = uigetfile(...
            {'*-data-??-???-????-??-??-??.mat','project data files'},...
            'select project data file');




        if fileName==0
            if GUI
                h = errordlg('no data loaded');
                uiwait(h)
                idlist = [];
                return
            else
                error('no data loaded')
            end
        end

    case 2 % assume we are (or are pointed) to the right directory
        f = searchFiles('.*-data-\d\d-.*\.mat','log',pathName);
        if isempty(f)
            if GUI
                h = errordlg('file not found');
                uiwait(h)
                idlist = [];
                return
            else
                error('file not found')
            end
        end
        fileName = f{1};



end

%=======================


%=================================================
% LOAD DATA FILE AND ASSIGN INDIVIDUAL VARIABLES
%=================================================

data = load([pathName,fileName]); %loads everything into the structure data


% -------- load idlist -------------
%find which idlists there are
dataFieldNames = fieldnames(data);
idnameListIdx = strmatch('idlist',dataFieldNames);
idnameList = dataFieldNames(idnameListIdx);

% check whether we want to forceNew
if forceNew
    % look for '2' in name, remove all that don't have any. If idlist
    % becomes empty, return idlist=-1.
    twoCell = regexp(idnameList,'2');
    noTwoIdx = cellfun(@isempty,twoCell);
    idnameList(noTwoIdx) = [];

    if isempty(idnameList)
        continueLoadIdlist = 0;
        idlist = -1;
    else
        continueLoadIdlist = 1;
    end
else
    continueLoadIdlist = 1;
end

if continueLoadIdlist

    idIdx = [];
    if ~isempty(idname)
        if strcmp(idname,'last')
            idIdx = find(strcmp(idnameList,data.lastResult));
            if ~isempty(idIdx)
                idname = idnameList{idIdx};
            end
        else
            idIdx = find(strcmp(idnameList,idname));
        end
    end

    if isempty(idIdx)
        %have the user choose, if there is more than one entry left
        switch length(idnameList)
            case 0 %no idlist loaded. if GUI, continue w/o loading

                if GUI
                    h = warndlg('no idlist found in data');
                    uiwait(h);
                    idname = '';
                else
                    error('no idlist found in data')
                end

            case 1 %only one idlist loaded. Continue

                idname = char(idnameList);

            otherwise %let the user choose
                idSelect = chooseFileGUI(idnameList);
                if isempty(idSelect)
                    idname = '';
                else
                    idname = idnameList{idSelect};
                end
        end
    end
    if isempty(idname)
        if GUI
            % continue loading if more than just idlist
            idlist = [];
            if nargout == 1
                return
            end
        else
            error('file not found')
        end
    else
        idlist = data.(idname);
        % store idname in idlist
        idlist(1).stats.idname = idname;
    end

end


%------------ dataProperties ------------------
if nargout > 1
    if ~isfield(data,'dataProperties')
        if GUI
            h = errordlg('No dataProperties in project data: corrupt data file');
            uiwait(h)
            return
        else
            error('No dataProperties in project data: corrupt data file');
        end
    else
        dataProperties = data.dataProperties;
    end
end


%------------ projectProperties
if nargout > 2
    if ~isfield(data,'projProperties')
        if GUI
            h = errordlg('No projProperties in project data!');
            uiwait(h)
            return
        else
            error('No projProperties in project data!');
        end
    else
        projectProperties = data.projProperties;
    end
end

% ---------- slist
if nargout > 3
    if ~isfield(data,'slist')
        if GUI
            h = errordlg('No slist in project data!');
            uiwait(h)
            return
        else
            error('No slist in project data!')
        end
    else
        slist = data.slist;
    end
end


%--------- movie

if nargout > 4

    switch noMovie
        case 0

            %--------------try to load filtered movie
            %try to find filenames in the path from which projectData has been loaded
            filteredMovie = cdLoadMovie('latest',pathName);

            %test if everything correctly loaded
            if ~exist('filteredMovie','var')
                if GUI
                    h = errordlg('no movie found');
                    uiwait(h)
                    return
                else
                    error('no movie found')
                end
            end

        case 1
            % loadStruct for filtered movie
            [dummy,dummy,loadStruct] = cdLoadMovie('filtered',pathName,-1);
            filteredMovie = loadStruct;
        case 2
            % loadStruct for corr/raw
            [dummy,dummy,loadStruct] = cdLoadMovie('corr/raw',pathName,-1);
            filteredMovie = loadStruct;
            
        case 3
            % return empty
            filteredMovie = [];

    end

end

