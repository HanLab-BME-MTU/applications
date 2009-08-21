function [projList,projPaths]=getProj(varargin)
% GETPROJ returns paths to roi_x and sub_x projects matching user query

% INPUT: One or more comma-separated search strings used to search all
%        roi_x and sub_x directory paths found in a user-selected parent
%        directory. The search is NOT case-sensitive.
%        Input can also be a cell array containing strings.
%
%        Note: if pwd is given as input, user is not asked for parent
%        directory; instead the working directory is used. In this case all
%        projects are returned. This option is not compatible with user
%        queries to narrow down projects.
%
% OUTPUT: projList       : structure containing fields:
%                             .anDir - path to analysis folder
%                             .imDir - path to image folder
%                          for each project fitting the query.
%                          this is saved in the top-level directory.
%         projPaths      : cell array containing roi/subroi paths only




projList=[];
projPaths=[];

if nargin>=1 && isempty(varargin{1,1})
    varargin=[];
end

n=[];
nStr=0;
if ~isempty(varargin)

    % check whether input is a cell array
    nStr=nargin;
    if isempty(varargin{1,1})
        nStr=0;
    end
    if iscell(varargin{1,1})
        varargin=varargin{1,1}';
        nStr=length(varargin);
    end

    % check whether inputs are strings
    inputStrings=cellfun(@(y) ischar(y), varargin);
    if sum(inputStrings)~=nStr
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
if nStr==1 && isequal(varargin{1},pwd)
    topDir=pwd;
else
    topDir=uigetdir(pwd,'Please select parent directory containing projects');
end
if topDir==0
    return
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
projPaths=[];

if isempty(roiDirList) && isempty(subDirList)
    h=msgbox('No projects found.');
    uiwait(h);
    return
end

% check roi list for matches to input strings
projCount=0;
if ~isempty(roiDirList)
    tempROI=ones(length(roiDirList),1);
    for i=1:nStr
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
    for i=1:nStr
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
    h=msgbox('No projects found.');
    uiwait(h);
    return
end

% convert structure to cell to view easier - but don't save it.
projPaths=struct2cell(projList);
projPaths=projPaths(2,:)';

save([topDir filesep 'projList' n],'projList')

