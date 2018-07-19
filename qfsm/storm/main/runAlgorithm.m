function runAlgorithm(useIncrementalAlgorithm)

% Open biggest possible matlabpool
if matlabpool('size') == 0 % No matlabpool open
    success = false;
    for nWorkers=[12,8,4]
        try
            matlabpool('open','local',nWorkers);
            success = true;
            break; 
        end
    end
    if success
        fprintf('Main: Opened matlabpool with %d workers!\n',matlabpool('size'))
    else
        disp('Main: Could not open matlabpool!')
    end
else
    fprintf('Main: Using existing matlabpool with %d workers!\n',matlabpool('size'))
end

% Define the root path
rootPath = [getStormPath() filesep '_queue' filesep];

% Find the first directory ending with '+'
list = dir(rootPath);
for i=3:numel(list)
    if strcmp(list(i).name(end),'+')
        dirName = list(i).name;
        break;
    end
end

% Run the algorithm if a directory has been found
if exist('dirName','var')
    dirPath = [rootPath dirName filesep];
    if nargin == 1
        disp('Main: Using incremental algorithm!');
        exitflag = algorithmIncremental(dirPath);
    else
        exitflag = algorithm(dirPath);
    end
    fprintf('%s\n',exitflag);
else
    disp('Main: No directory found!');
end

end

