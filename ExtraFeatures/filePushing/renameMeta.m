function [ output_args ] = renameMeta(groupList,beforeRename,afterRename)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(groupList)
    
    fileNameBeforeRename = [char(groupList(i,2)) filesep beforeRename] ;
    fileNameAfterRename = [char(groupList(i,2)) filesep afterRename];
    
    movefile(fileNameBeforeRename, fileNameAfterRename);
    

end

