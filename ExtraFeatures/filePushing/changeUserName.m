function [ groupList] = changeUserName(saveDir,filename,name2replace,userName,groupList)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if isempty(saveDir);
saveDir = uigetdir(pwd,'Select output directory for renamed groupList');
end 

    
 
for iGroup = 1:length(groupList)
 x  = regexprep(char(groupList(iGroup,2)),name2replace,userName);
 groupList{iGroup,2} = x;
end

fullFilename = [saveDir filesep filename];
save(fullFilename,'groupList');
end
