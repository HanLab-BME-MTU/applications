function [ projList] = changeUserNameProjList(saveDir,filename,name2replace,userName,projList)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if isempty(saveDir);
saveDir = uigetdir(pwd,'Select output directory for renamed projList');
end 

    
 
for iProj = 1:length(projList)
 projList(iProj).imDir = regexprep(projList(iProj).imDir,name2replace,userName);
 projList(iProj).anDir = regexprep(projList(iProj).anDir,name2replace,userName);
 
end

fullFilename = [saveDir filesep filename];
save(fullFilename,'projList');
end
