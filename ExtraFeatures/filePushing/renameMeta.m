function [ output_args ] = renameMeta(projList,beforeRename,afterRename)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for i = 1:length(projList)
    
    fileNameBeforeRename = [projList(i).anDir filesep beforeRename] ;
    fileNameAfterRename = [projList(i).anDir filesep afterRename];
    
    movefile(fileNameBeforeRename, fileNameAfterRename);
    

end

