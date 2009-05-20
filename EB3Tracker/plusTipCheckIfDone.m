function [result,notDone]=plusTipCheckIfDone
% plusTipCheckIfDone gives projData timestamp for each project and list of unanalyzed projects
%
% synopsis: [result,notDone]=plusTipCheckIfDone
%
result=[];
notDone=[];

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
allProjList=sort(allProjList);
allProjList=cellfun(@(i) formatPath(i),allProjList,'uniformoutput',0);

% list of all the files in each project directory
dirContents=cellfun(@(i) dir([i filesep 'meta']),allProjList,'uniformoutput',0);
% convert each structure to a cell array
dirContentsCell=cellfun(@(i) struct2cell(i),dirContents,'uniformoutput',0);
% get index for projData file, if it exists
projListLoc=cellfun(@(i) find(strcmpi(i(1,:),'projData.mat')),dirContentsCell,'uniformoutput',0);
% get timestamp for each projList file
dates=cellfun(@(i,j) i(2,j),dirContentsCell,projListLoc,'uniformoutput',0); 

% find out if any of them don't have dates
notDone=find(cell2mat(cellfun(@(i) isempty(i),dates,'uniformoutput',0)));

result=[allProjList dates];