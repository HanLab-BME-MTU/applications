function [projectDirectories]=projList2Mat(projList)
% return list of paths for project analysis directories in projList struct

if isempty(projList)
    projectDirectories=[];
else
    projectDirectories=struct2cell(projList);
    if isempty(strfind(projectDirectories{2,1},[filesep 'images']))
        projectDirectories=projectDirectories(2,:)';
    else
        projectDirectories=projectDirectories(1,:)';
    end
    homeDir=pwd;
    projectDirectories=cellfun(@(i) formatPath(i),projectDirectories,'uniformoutput',0);
    cd(homeDir)
end