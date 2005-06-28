function [idlist,dataProperties,projectProperties,slist,filteredMovie] = loadProjectData(fileName,pathName,idname)
%LOADPROJECTDATA loads experimental Chromdyn data from file
%
% SYNOPSIS [idlist,dataProperties,projectProperties,slist,filteredMovie] = loadProjectData(fileName,pathName,idname)
%
% INPUT    fileName  (opt) name of data file. if 1, the program will load
%                       GUI mode
%          pathName  (opt) name of path for data file. If fileName, but no
%                       pathName is specified, currentDir is used
%          idname    (opt) name of idlist to be loaded. Can be "last"
%
% OUTPUT   idlist           user selected idlist if several possible
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

if nargin == 0 || isempty(fileName)
    loadData = 1;
    GUI = 0;
else
    if fileName == 1
        loadData = 1;
        GUI = 1;
    else
        loadData = 0;
        GUI = 0;
    end
end

if nargin < 2 || isempty(pathName)
    pathName = pwd;
end
if ~strcmp(pathName(end),filesep)
    pathName = [pathName,filesep];
end

if nargin < 3 || isempty(idname)
    idname = [];
end

%=================


%========================
% FIND FILE IF NECESSARY
%========================

if loadData
    oldDir = pwd;
    
    %check if default biodata-dir exists and cd if exist
    mainDir = cdBiodata(2);
    
    %get project data file
    [fileName,pathName] = uigetfile({'*-data-??-???-????-??-??-??.mat','project data files'},'select project data file');
    
    
    
    
    if fileName==0
        if GUI
            h = errordlg('no data loaded')
            uiwait(h)
            return
        else
            error('no data loaded')
        end
    end
    
else
    
    if ~exist(fileName) && ~ exist([pathName,fileName])
        if GUI
            h = errordlg('file not found')
            uiwait(h)
            return
        else
            error('file not found')
        end
    end
    
end

%=======================


%=================================================
% LOAD DATA FILE AND ASSIGN INDIVIDUAL VARIABLES
%=================================================

data = load([pathName,filesep,fileName]); %loads everything into the structure data


% -------- load idlist -------------
%find which idlists there are
dataFieldNames = fieldnames(data);
idnameListIdx = strmatch('idlist',dataFieldNames);
idnameList = dataFieldNames(idnameListIdx);

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
                h = warndlg('no idlist found in data')
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
        % continue loading
        idlist = [];
    else
        error('file not found')
    end
else
    idlist = data.(idname);
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
            h = errordlg('No slist in project data!')
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
    
    %--------------try to load filtered movie
    %try to find filenames in the path from which projectData has been loaded
    filteredMovie = cdLoadMovie('latest',pathName);
    
    %test if everything correctly loaded
    if ~exist('filteredMovie','var')
        if GUI
            h = errordlg('file not found')
            uiwait(h)
            return
        else
            error('file not found')
        end
        error('no movie found')
        return
    end
    
end

