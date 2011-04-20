function [ output_args ] = renameMeta(groupList)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(groupList)
    
    fileNameBeforeRename = [char(groupList(i,2)) filesep 'meta'] ;
    fileNameAfterRename = [char(groupList(i,2)) filesep 'metaOld'];
    
    movefile(fileNameBeforeRename, fileNameAfterRename);
    

end

