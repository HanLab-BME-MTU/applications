function computeTMDelay(varargin)

if nargin > 0
    rootDirectory = varargin{1};
else
    % Ask for the root directory.
    rootDirectory = uigetdir('', 'Select a root directory:');
end

if ~ischar(rootDirectory)
    return;
end

if nargin == 2
    redo = varargin{2};
else
    redo = 1;
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
    movieName = ['Movie ' num2str(iMovie) '/' num2str(numel(paths))];    
    path = paths{iMovie};
    
    % Check is the computation has been already performed.    
    tmDelayFileName = [path filesep 'tmDelay.mat'];
    
    if ~redo && exist(tmDelayFileName, 'file')
        disp([movieName ': Analysis has be already done. (SKIPPING)']);
        continue;
    end
    
    % Get the 2 FSMCenter directory paths (Actin, TM2 or TM4)
    subDirNames = {'Actin', 'TM2', 'TM4'};
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
    
    % Get the mask filename from the first FSMCenter directory (Actin
    % directory that is).
    maskFileNames = dir([subDirPaths{1} filesep 'edge' filesep...
        'cell_mask' filesep 'mask*.tif']);
    numFrames = numel(maskFileNames);
    if ~numFrames
        disp([movieName ': Unable to locate mask files. (SKIPPING)']);
        continue;
    end
    
    % Compute the distance transforms
    D = cell(numel(maskFileNames), 1);    
    h = waitbar(0, [movieName ': Compute distance transforms...']);
    for iFrame = 1:numFrames
        fileName = [subDirPaths{1} filesep 'edge' filesep 'cell_mask'...
            filesep maskFileNames(iFrame).name];
        
        BW = imread(fileName);

        % Compute the distance transform
        D{iFrame} = bwdist(max(BW(:)) - BW);
        
        waitbar(iFrame / numel(maskFileNames), h);        
    end
    
    % Load MPM files
    MPMs = cell(2, 1);
    
    for k = 1:2
        fileName = [subDirPaths{k} filesep 'tack' filesep 'mpm.mat'];
    
        if ~exist(fileName, 'file')
            disp([movieName ': Unable to find the mpm.mat file in '...
                subDirPaths{k} '. (SKIPPING)']);
            continue;
        end
    
        load(fileName);
    
        if ~exist('MPM', 'var')
            disp([movieName ': Unable to find the Actin MPM variable in '...
                subDirPaths{k} '/tack/mpm.mat. (SKIPPING)']);
            continue;
        end
        
        MPMs{k} = MPM;
        clear M MPM;
    end
    
    if size(MPMs{1}, 2) ~= size(MPMs{2}, 2) && size(MPMs{1}, 2) ~= 0
        disp([movieName ': Number of frames in MPM variables differ '...
            '(SKIPPING MOVIE)']);
        continue;
    end
    
    % Compute the delay between TM birth and Actin birth
    
    tmNumBirths = 0;
    % tmBirths = [y, x, distToEdge, meanDelayFromActinBirths];
    tmBirths = zeros(numel(find(MPMs{2} == 0)) / 2 + 1, 6);
    waitbar(0, h, [movieName ': Compute TM delays...']);
    
    % First frame:
    waitbar(1 / numFrames, h);
    for i = 1:size(MPMs{2}, 1)
        if MPMs{2}(i, 1)
            tmNumBirths = tmNumBirths + 1;            
            p = MPMs{2}(i, 1:2);            
            tmBirths(tmNumBirths, 1:2) = p;
            tmBirths(tmNumBirths, 3) = D{1}(p(1), p(2));
            
            % Find actin tracks for which distance to p <= 4
            distToActinTracks = createDistanceMatrix(p, MPMs{1}(:, 1:2));
            iClosestActinTracks = find(distToActinTracks <= 4);
            if numel(iClosestActinTracks)
                tmBirths(tmNumBirths, 4) = 0;
            else
                tmBirths(tmNumBirths, 4) = -1;
            end
        end
    end
    % Rest of frames:
    for j = 3:2:size(MPMs{2}, 2)
        t = ceil(j / 2);
        waitbar(t / numFrames, h);
        for i = 1:size(MPMs{2}, 1)
            if ~MPMs{2}(i, j - 1) && MPMs{2}(i, j)
                tmNumBirths = tmNumBirths + 1;
                p = MPMs{2}(i, j:j+1);
                tmBirths(tmNumBirths, 1:2) = p;
                tmBirths(tmNumBirths, 3) = D{t}(p(1), p(2));

                % Find the closest Actin track to p at frame ceil(j/2)
                distToActinTracks = createDistanceMatrix(p, MPMs{1}(:, j:j+1));
                iClosestActinTracks = find(distToActinTracks <= 4);
                if numel(iClosestActinTracks)
                    meanDelay = 0;
                    for k = 1:numel(iClosestActinTracks)
                        j2 = j;
                        while j2 > 2 && MPMs{1}(iClosestActinTracks(k), j2 - 2)
                            j2 = j2 - 2;
                        end
                        meanDelay = meanDelay + t - ceil(j2/2);
                    end
                    tmBirths(tmNumBirths, 4) = meanDelay / numel(iClosestActinTracks);
                else
                    tmBirths(tmNumBirths, 4) = -1;
                end
            end
        end
    end
    close(h);
    tmBirths = tmBirths(1:tmNumBirths, :); %#ok<NASGU>

    % Save delays.
    save(tmDelayFileName, 'tmBirths');
    clear tmBirths;
    
    % Clear all variables
    clear D;
    
    disp([movieName ': DONE']);
end

toc
end