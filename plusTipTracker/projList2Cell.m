function [projectDirectories]=projList2Cell(projList)
% return nx2 cell array of [project images] paths from projList struct

if isempty(projList)
    projectDirectories=[];
else   
    projectDirectories=struct2cell(projList)';
    if isempty(strfind(formatPath(projectDirectories{1,2}),[filesep 'images']))
        projectDirectories=projectDirectories(:,2:-1:1);
    end
end
