function [groupList]=combineGroupListFiles(saveResult)
% combineGroupListFiles allows selection of multiple groupList files to combine
%
% SYNOPSIS: [groupList]=combineGroupListFiles(saveResult)
%
% INPUT   : saveResult: 1 to save new groupList file, 0 to return it to
%                       workspace only
% OUTPUT  : groupList  : see plusTipPickGroups.m for details
%
% 2009/08 - Kathryn Applegate, Matlab R2008a

if nargin<1
    saveResult=0;
end

groupList=[];
temp=[];
userEntry='Yes';
while strcmp(userEntry,'Yes')
    [fileName,pathName] = uigetfile('*.mat','Select groupList.mat file');
    if fileName==0
        return
    end
    load([pathName filesep fileName]);
    temp=[temp; groupList];
    disp(['Selected: ' pathName fileName])
    userEntry = questdlg('Select another groupList.mat file?');
end
clear groupList;
groupList=temp;

% format path for current OS
nProj=length(groupList);
curDir=pwd;
groupList(:,2)=cellfun(@(x) formatPath(groupList{x,2}),mat2cell([1:nProj]',ones(nProj,1),1),'uniformOutput',0);
cd(curDir)

if saveResult==1
    temp=inputdlg({'Enter file name:'},'',1);
    dirName=uigetdir(pwd,'Select output directory.');
    save([dirName filesep temp{1}],'groupList')
end