function [projList,projPathParsed]=getProj(varargin)
% GETPROJ returns paths to roi_x and sub_x projects matching user query

% INPUT: one or more comma-separated search strings used to search all
%        roi_x and sub_x directory paths found in a user-selected top-level
%        directory. the search is NOT case-sensitive.
%
%        Note: if pwd is given as input, user is not asked for top-level
%        directory; instead the working directory is used. in this case all
%        projects are returned. this option is not compatible with user
%        queries to narrow down projects.
%
% OUTPUT: projList       : structure containing fields:
%                             .anDir - path to analysis folder
%                             .imDir - path to image folder
%                          for each project fitting the query.
%                          this is saved in the top-level directory.
%         projPathParsed : cell array containing easier-to-read info about
%                          the projects


% check whether inputs are strings
n=[];
if ~isempty(varargin)
    inputStrings=cellfun(@(y) ischar(y), varargin);
    if sum(inputStrings)~=nargin
        error('input arguments must be strings')
    end

    % create string for naming projList with query
    for i=1:length(varargin)
        if strcmp(varargin{i},pwd)~=1
            n=[n '_' varargin{i}];
        end
    end
end

% if input is current working directory, look here and in all
% sub-directories for roi_x directories. otherwise, ask the user to select
% the top-level directory
if nargin==1 && isequal(varargin{1},pwd)
    topDir=pwd;
else
    topDir=uigetdir(pwd,'Please select top-level directory containing targets');
end


disp(topDir);
% get semi-colon- (windows) or colon- (linux) separated list of all
% sub-directories under the top directory
p=genpath(topDir);
if ispc
    tempDirList=strrep(p,';',' ');
else
    tempDirList=strrep(p,':',' ');
end
% get cell array of "roi_x" directories
roiDirList=regexp(tempDirList,['\S*roi_\d\s'],'match')';
subDirList=regexp(tempDirList,['\S*sub_\d\s'],'match')';

% initialize as empty in case no projects found
projList=[];
projPathParsed=[];

if isempty(roiDirList) && isempty(subDirList)
    disp('no projects found')
    return
end

% check roi list for matches to input strings
projCount=0;
if ~isempty(roiDirList)
    tempROI=ones(length(roiDirList),1);
    for i=1:nargin
        testStr = varargin{i};
        tempROI=tempROI & cellfun(@(y) ~isempty(strfind(y,lower(testStr))),lower(roiDirList));
    end

    matches=find(tempROI);

    for i=1:length(matches)
        roiDir=roiDirList{matches(i),1};
        projList(i,1).imDir=[roiDir(1:end-6) 'images'];
        projList(i,1).anDir=roiDir(1:end-1);
        projCount=projCount+1;
    end
end
% do this also for sub list
if ~isempty(subDirList)
    tempSUB=ones(length(subDirList),1);
    for i=1:nargin
        testStr = varargin{i};
        tempSUB=tempSUB & cellfun(@(y) ~isempty(strfind(y,lower(testStr))),lower(subDirList));
    end

    matches=find(tempSUB);

    for j=1:length(matches)
        subDir=subDirList{matches(j),1};
        homeDir=pwd;
        projList(projCount+j,1).anDir=subDir(1:end-1);
        cd(projList(projCount+j,1).anDir)
        cd(['..' filesep '..' filesep '..' filesep 'images'])
        projList(projCount+j,1).imDir=pwd;
        cd(homeDir)
    end
end

% return if nothing
if isempty(projList)
    disp('no projects found')
    return
end

projPathParsed{size(projList,1),6}=[];
for iProj=1:size(projList,1)

    currentROI=formatPath(projList(iProj,1).anDir);
    % parse the path to get "words" used to identify target, oligo,
    % movie, and roi
    nChar=length(currentROI);
    if ispc
        filesepLoc=regexp(currentROI,'\\');
    else
        filesepLoc=regexp(currentROI,'\/');
    end
    wordStart=[1 filesepLoc+1]; wordEnd=[filesepLoc-1 nChar];
    words=cell(length(wordStart),1);
    for iWord=1:length(wordStart)
        words{iWord,1}=currentROI(wordStart(iWord):wordEnd(iWord));
    end

    % assign values to cell
    if strfind(words{end},'roi')
        projPathParsed{iProj,1}=words{end-4,1}; % date
        projPathParsed{iProj,2}=words{end-3,1}; % target
        projPathParsed{iProj,3}=words{end-2,1}; % oligo
        projPathParsed{iProj,4}=words{end-1,1}; % movie number
        projPathParsed{iProj,5}=words{end-0,1}; % roi number
    else
        projPathParsed{iProj,1}=words{end-6,1}; % date
        projPathParsed{iProj,2}=words{end-5,1}; % target
        projPathParsed{iProj,3}=words{end-4,1}; % oligo
        projPathParsed{iProj,4}=words{end-3,1}; % movie number
        projPathParsed{iProj,5}=words{end-2,1}; % roi number
        projPathParsed{iProj,6}=words{end-0,1}; % sub-roi number
    end
end


save([topDir filesep 'projList' n],'projList')

