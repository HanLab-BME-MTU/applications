function makeImageDirectory
% MAKEIMAGEDIRECTORY moves tifs from parent to new subdirectory 'images'

% get parent directory from user containing multiple projects for which
% the action should be done
topDir=uigetdir(pwd,'Please select parent directory containing projects');

% look for tifs that do not include "overlay" in their name (overlay images
% are produced by kathryn's EB3 comet detection)
[listOfImages]=searchFiles('.tif','overlay',topDir,1);

% find the unique directory names
projDirList = unique(listOfImages(:,2));

% if the parent path includes "images," the operation won't be done; if
% not, a new directory called images is created and tifs are moved there
for i=1:length(projDirList)
    if isempty(findstr(projDirList{i}, 'images'))
        mkdir([projDirList{i} filesep 'images'])
        movefile([projDirList{i} filesep '*.TIF'],[projDirList{i} filesep 'images'])
    end
end