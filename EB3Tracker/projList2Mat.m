function [projectDirectories]=projList2Mat(projList)
% return list of paths for project analysis directories in projList struct

projectDirectories=struct2cell(projList);
projectDirectories=projectDirectories(2,:)';
homeDir=pwd;
projectDirectories=cellfun(@(i) formatPath(i),projectDirectories,'uniformoutput',0);
cd(homeDir)