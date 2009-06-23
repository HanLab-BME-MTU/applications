function [result,notDone]=plusTipCheckIfDone(fileType)
% plusTipCheckIfDone returns directory name and timestamp for creation of either
% movieInfo.mat (if fileType=1) or projData.mat (if fileType=2) for each
% project, as well as which projects have not been analyzed.

result=[];
notDone=[];

if nargin < 1
    disp('--checkIfDone: Incorrect number of input arguments!');
    errFlag  = 1;
    return
end

[FileName,PathName] = uigetfile('*.mat','Select projList.mat file');

% load projList
if FileName==0
    return
else
    homeDir=pwd;
    cd(PathName)
    projList=load([PathName filesep FileName]);
    projList=projList.projList;
    cd(homeDir)
end

% get sorted list of all projects
allProjList=struct2cell(projList);
allProjList=allProjList(2,:)';
%allProjList=sort(allProjList);
allProjList=cellfun(@(i) formatPath(i),allProjList,'uniformoutput',0);

switch fileType
    case 1

        % list of all the files in each project directory
        dirContents=cellfun(@(i) dir([i filesep 'feat']),allProjList,'uniformoutput',0);
        % convert each structure to a cell array
        dirContentsCell=cellfun(@(i) struct2cell(i),dirContents,'uniformoutput',0);
        % get index for movieInfo file, if it exists
        movieInfoLoc=cellfun(@(i) find(strcmpi(i(1,:),'movieInfo.mat')),dirContentsCell,'uniformoutput',0);
        % get timestamp for each movieInfo file
        dates=cellfun(@(i,j) i(2,j),dirContentsCell,movieInfoLoc,'uniformoutput',0);
    case 2
        % list of all the files in each project directory
        dirContents=cellfun(@(i) dir([i filesep 'feat']),allProjList,'uniformoutput',0);
        % convert each structure to a cell array
        dirContentsCell=cellfun(@(i) struct2cell(i),dirContents,'uniformoutput',0);
        % get index for projData file, if it exists
        projListLoc=cellfun(@(i) find(strcmpi(i(1,:),'projData.mat')),dirContentsCell,'uniformoutput',0);
        % get timestamp for each projList file
        dates=cellfun(@(i,j) i(2,j),dirContentsCell,projListLoc,'uniformoutput',0);
end

% find out if any of them don't have dates
notDone=find(cell2mat(cellfun(@(i) isempty(i),dates,'uniformoutput',0)));

result=[allProjList dates];