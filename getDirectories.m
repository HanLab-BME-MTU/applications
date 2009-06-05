function paths = getDirectories(rootDirectory)
% This function finds every path which contains an Actin and TM*
% subdirectories from rootDirectory.

% Get the list of all subdirectories of rootDirectory (including itself).

if ~exist(rootDirectory, 'dir')
    error([rootDirectory 'is not a valid directory.']);
end

paths = {rootDirectory};
iPath = 1;

while iPath <= numel(paths)
    contents = dir(paths{iPath});

    for i=1:numel(contents)
        name = contents(i).name;
    
        if contents(i).isdir &&...
                ~(strcmp(name, '.') || strcmp(name, '..'))
            containsActin = ~isempty(regexpi(name, 'actin'));
            containsTM = ~isempty(regexpi(name, 'TM'));
            
            if containsActin == containsTM
                paths = cat(1, paths, [paths{iPath} filesep name]);
            end
        end
    end
    
    iPath = iPath + 1;
end

% Select paths which contains an Actin and TM* subdirectories

selectedPaths = {};

for iPath = 1:numel(paths)
    contents = dir(paths{iPath});
    
    actinFound = 0;
    TMFound = 0;
    for i = 1:numel(contents)
        name = contents(i).name;
        if contents(i).isdir
            containsActin = ~isempty(regexpi(name, 'actin'));
            containsTM = ~isempty(regexpi(name, 'TM'));
            
            % to cope with directory name containing both 'actin' and 'TM'
            % strings.
            if xor(containsActin, containsTM)        
                actinFound = actinFound | containsActin;
                TMFound = TMFound | containsTM;
            end
        end
    end
    if actinFound && TMFound
        selectedPaths = cat(1, selectedPaths, paths{iPath});
    end
end

paths = selectedPaths;


end