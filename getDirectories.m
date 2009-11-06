function paths = getDirectories(rootDirectory)
% This function finds every path which contains either Actin and TMx
% subdirectories or TMx TMy subdirectories from rootDirectory. These 2
% subdirectories must be valid FSMCenter project directories.

if ~exist(rootDirectory, 'dir')
    error([rootDirectory ' is not a valid directory.']);
end

disp('Find TMx/TMy/Actin subdirectories. Please wait...');
str = genpath(rootDirectory);
paths = textscan(str, '%s', 'delimiter', pathsep);
paths = paths{1};
finalPaths = {};

for iPath = 1:numel(paths)
    contents = dir(paths{iPath});
    
    found = 0;
    
    for i = 1:numel(contents)
        name = upper(contents(i).name);
        
        if contents(i).isdir &&  ...
                ~(strcmp(name, '.') || strcmp(name, '..')) && ...
                ~isempty(strmatch(name, {'ACTIN', 'TM2', 'TM4', 'TM5NM1'}, 'exact'))
            % TODO: Check that name is a valid FSMCEnter project directory
            found = found + 1;
        end
    end
    
    if found == 2
        finalPaths = cat(1, finalPaths, paths{iPath});
    end
end

paths = finalPaths;
