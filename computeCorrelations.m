function computeCorrelations(varargin)
% Compute correlation of map ch1 and map ch2

if nargin >= 1
    rootDirectory = varargin{1};
else
    % Ask for the root directory.
    rootDirectory = uigetdir('', 'Select a root directory:');

    if ~ischar(rootDirectory)
        return;
    end
end

if nargin >= 2
    redo = varargin{2};
else
    redo = 0;
end

% Get all subdirectories containing Actin & TM.
paths = getDirectories(rootDirectory);

disp('List of directories:');

for iMovie = 1:numel(paths)
    disp([num2str(iMovie) ': ' paths{iMovie}]);
end

disp('Process all directories...');

tic

for iMovie = 1:numel(paths)
    movieName = ['Movie ' num2str(iMovie) '/'...
        num2str(numel(paths))];
    
    path = paths{iMovie};
    
    % Check whether the computation has already been done.
    outputDirName = [path filesep 'correlations'];
    
    if exist(outputDirName, 'dir')
        if ~redo
            disp([movieName ': DONE']);
            continue;
        end
    else
        mkdir(outputDirName);
    end
    
    % Get the 2 FSMCenter directory paths (Actin, TM2 or TM4 or TM5NM1)
    % Linux is case sensitive on file system so gives every posibility
    if isunix
        subDirNames = {'actin', 'Actin', 'TM2', 'TM4', 'TM5NM1'};
    else
        subDirNames = {'actin', 'TM2', 'TM4', 'TM5NM1'};
    end
    subDirPaths = cell(2, 1);
    indSubDir = 0;
    for i = 1:numel(subDirNames)
        subDirPath = [path filesep subDirNames{i}];
        if exist(subDirPath, 'dir')
            indSubDir = indSubDir + 1;
            subDirPaths{indSubDir} = subDirPath;
        end
    end
    
    if indSubDir ~= 2
        disp([movieName ': Unable to find the 2 FSM Center directories. (SKIPPING)']);
        continue;
    end
    
    speedMapFileNames = cell(2, 1);
    error = 0;
    
     for k = 1:2
        speedMapFileNames{k} = dir([subDirPaths{k} filesep 'post' ...
            filesep 'mat' filesep 'speedMap*.mat']);
        
        if isempty(speedMapFileNames{k})
            disp([movieName ': Unable to locate speed map files. (SKIPPING)']);
            error = 1;
            break;
        end
     end
     
     if error
         continue;
     end
     
     for iFile = 1:numel(speedMapFileNames{1})
         fileName = [subDirPaths{1} filesep 'post' filesep 'mat'...
             filesep speedMapFileNames{1}(iFile).name];
         [dummy1, dummy2, no] = getFilenameBody(fileName);
         
         load(fileName);
         speedMap = speedMap / max(speedMap(:));
         speedMap1 = speedMap;
         
          fileName = [subDirPaths{2} filesep 'post' filesep 'mat'...
             filesep speedMapFileNames{2}(iFile).name];
         load(fileName);
         speedMap = speedMap / max(speedMap(:));
         speedMap2 = speedMap;
         
         speedCorr = (speedMap1 - speedMap2) / sqrt(2); %#ok<NASGU>
         save([outputDirName filesep 'speedCorr_' no '.mat'], 'speedCorr');
     end
end

toc;

end

