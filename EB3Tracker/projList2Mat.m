function [projectDirectories]=projList2Mat(projList)
% return list of paths for project analysis directories in projList struct

projectDirectories=struct2cell(projList);
projectDirectories=projectDirectories(2,:)';
projectDirectories=cellfun(@(i) formatPath(i),projectDirectories,'uniformoutput',0);