function [projList]=combineProjListFiles(saveResult)
% combineProjListFiles allows selection of multiple projList files to combine

% SYNOPSIS: [projList]=combineProjListFiles(saveResult)
%
% INPUT   : saveResult: 1 to save new projList file, 0 to return it to
%                       workspace only
% OUTPUT  : projList  : see getProj.m for details
%
% 2009/08 - Kathryn Applegate, Matlab R2008a

if nargin<1
    saveResult=0;
end
projList=[];
temp=[];
userEntry='Yes';
while strcmp(userEntry,'Yes')
    [fileName,pathName] = uigetfile('*.mat','Select projList.mat file');
    if fileName==0
        return
    end
    load([pathName filesep fileName]);
    temp=[temp; projList];
    disp(['Selected: ' pathName filesep fileName])
    userEntry = questdlg('Select another projList.mat file?');
end
clear projList;
projList=temp;

if saveResult==1
    temp=inputdlg({'Enter file name:'},'',1);
    dirName=uigetdir('Select output directory.');
    save([dirName filesep temp{1}],'projList')
end