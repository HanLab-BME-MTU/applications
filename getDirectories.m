function paths = getDirectories(rootDirectory, numSubDirs, subDirList, varargin)
% This function finds every path which contains 'numbSubDirs' sub folders
% among the list of possible names subDirList. A function handler can be
% specified to validate each sub directory.

if ~exist(rootDirectory, 'dir')
    error([rootDirectory ' is not a valid directory.']);
end

if numSubDirs <= 0 || numSubDirs > numel(subDirList)
    error('invalid numSubDirs parameter.');
end

subDirList = upper(subDirList);

if nargin < 4 || isempty(varargin{1})
    isValidSubDir = @(x) true;
else
    isValidSubDir = varargin{1};
end

disp('Find sub directories. Please wait...');
str = genpath(rootDirectory);
paths = textscan(str, '%s', 'delimiter', pathsep);
paths = paths{1};
finalPaths = {};

for iPath = 1:numel(paths)
    contents = dir(paths{iPath});
    
    found = 0;
    
    for i = 1:numel(contents)
        name = contents(i).name;
        
        if contents(i).isdir &&  ...
                ~(strcmp(name, '.') || strcmp(name, '..')) && ...
                ~isempty(strmatch(upper(name), subDirList, 'exact')) && ...
                isValidSubDir([paths{iPath} filesep name ])
            found = found + 1;
        end
    end
    
    if found == numSubDirs
        finalPaths = cat(1, finalPaths, paths{iPath});
    end
end

paths = finalPaths;
