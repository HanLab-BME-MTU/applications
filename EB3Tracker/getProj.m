function [projList]=getProj(varargin)
% GETPROJ returns paths to projects matching user query

% INPUT: one or more comma-separated search strings used to search all roi_x
% directories in a user-selected top level directory
%
% OUTPUT:
% projList: structure containing image and analysis directories for each
% sub-project fitting the query. this is saved in the top-level directory

if ~isempty(varargin)
    inputStrings=cellfun(@(y) ischar(y), varargin);
    if sum(inputStrings)~=nargin
        error('input arguments must be strings')
    end
end

topDir=uigetdir(pwd,'Please select top-level directory containing targets');
p=genpath(topDir);
tempDirList=strrep(p,';',' ');
roiDirList = regexp(tempDirList,'\S*\\roi_\d\s','match')'; % cell array of "roi_x" directories

temp=ones(length(roiDirList),1);
for i=1:nargin
    testStr = varargin{i};
    temp=temp & cellfun(@(y) ~isempty(strfind(y,testStr)),lower(roiDirList));
end

matches=find(temp);
for i=1:length(matches)
    roiDir=roiDirList{matches(i),1};
    projList(i,1).imDir=[roiDir(1:end-6) 'images'];
    projList(i,1).anDir=roiDir(1:end-1);
end

save([topDir filesep 'projList'],'projList')

