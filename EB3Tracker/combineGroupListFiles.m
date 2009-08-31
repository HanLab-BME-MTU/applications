function [projGroupDir,projGroupName]=combineGroupListFiles(saveResult)
% combineGroupListFiles allows selection of multiple groupList files to combine

% SYNOPSIS: [projGroupDir,projGroupName]=combineGroupListFiles(saveResult)
%
% INPUT   : saveResult: 1 to save new groupList file, 0 to return it to
%                       workspace only
% OUTPUT  : projGroupDir,projGroupName: see plusTipPickGroups.m for details
%
% 2009/08 - Kathryn Applegate, Matlab R2008a

if nargin<1 || isempty(saveResult)
    saveResult=0;
else
    saveDir=uigetdir(pwd,'Select output directory.');
end

homeDir=pwd;

groupList=[];
tempDir=[];
tempName=[];
h=msgbox('Select groupList files in order');
uiwait(h)
userEntry='Yes';
count=1;
while strcmp(userEntry,'Yes')
    [fileName,pathName] = uigetfile('*.mat','Select groupList.mat file');
    if fileName==0
        return
    end
    cd(pathName)
    load([pathName filesep fileName]);
    tempDir=[tempDir; projGroupDir];
    tempName=[tempName; projGroupName];
    tempFileName{count}=fileName;
    disp(['Selected: ' pathName filesep fileName])
    userEntry = questdlg('Select another groupList.mat file?');
    count=count+1;
end

clear projGroupDir projGroupName;
projGroupDir=tempDir;
projGroupName=tempName;

if saveResult==1
    fileName='groupList';
    [nameList,m,n]=unique(projGroupName);
    [b,idx]=sort(m);
    nameList=nameList(idx);
    for i=1:length(nameList)
        fileName=[fileName '_' nameList{i}];
    end
    save([saveDir filesep fileName],'projGroupName','projGroupDir')
end
cd(homeDir)