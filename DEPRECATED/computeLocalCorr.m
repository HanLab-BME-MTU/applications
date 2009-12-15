function computeLocalCorr(varargin)
%      COMPUTELOCALCORR(forceRedo) Compute the correlation
%      between average speed and kinScore, along TM and Actin tracks.
%      'forceRedo' forces to recompute results.

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
    outputFileName = [path filesep 'localCorrelations.mat'];
    
    if exist(outputFileName, 'file') && ~redo
        disp([movieName ': DONE']);
        continue;        
    end
    
    % Get the 2 FSMCenter directory paths (Actin, TM2 or TM4)
    % Linux is case sensitive on file system so gives every posibility
    if isunix
        subDirNames = {'actin', 'Actin', 'TM2', 'TM4'};
    else
        subDirNames = {'actin', 'TM2', 'TM4'};
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
    
    % Get the mask filename from the first FSMCenter directory (Actin
    % directory that is).
    maskFileNames = dir([subDirPaths{1} filesep 'edge' filesep...
        'cell_mask' filesep 'mask*.tif']);
    numFrames = numel(maskFileNames);
    if ~numFrames
        disp([movieName ': Unable to locate mask files. (SKIPPING)']);
        continue;
    end
    
    % Get the list of speed, poly and depoly map file names.
    speedMapFileNames = cell(2, 1);
    polyMapFileNames = cell(2, 1);
    depolyMapFileNames = cell(2, 1);
    
    error = 0;
    
    for k = 1:2
        speedMapFileNames{k} = dir([subDirPaths{k} filesep...
            'post' filesep 'mat' filesep 'speedMap*.mat']);
        
        if isempty(speedMapFileNames{k})
            disp([movieName ': Unable to locate speed map files. (SKIPPING)']);
            error = 1;
            break;
        end
        
        polyMapFileNames{k} = dir([subDirPaths{k} filesep...
            'post' filesep 'mat' filesep 'polyMap*.mat']);
        
        if isempty(polyMapFileNames{k})
            disp([movieName ': Unable to locate poly map files. (SKIPPING)']);
            error = 1;
            break;
        end
        
        depolyMapFileNames{k} = dir([subDirPaths{k} filesep...
            'post' filesep 'mat' filesep 'depolyMap*.mat']);
        
        if isempty(polyMapFileNames{k})
            disp([movieName ': Unable to locate depoly map files. (SKIPPING)']);
            error = 1;
            break;
        end
    end
    
    if error
        continue;
    end

    % Load TM MPM file
    fileName = [subDirPaths{2} filesep 'tack' filesep 'mpm.mat'];
        
    if ~exist(fileName, 'file')
        disp([movieName ': Unable to find the mpm.mat file in '...
            subDirPaths{2} '. (SKIPPING)']);
        continue;
    end
    
    load(fileName);
    
    if ~exist('MPM', 'var')
        disp([movieName ': Unable to find the Actin MPM variable in '...
            subDirPaths{2} '/tack/mpm.mat. (SKIPPING)']);
        continue;
    end
    
    clear M;
    
    % Compute the distance transforms
    D = cell(numel(maskFileNames), 1);
    h = waitbar(0, [movieName ': Compute distance transforms...']);
    for iFrame = 1:numFrames
        fileName = [subDirPaths{1} filesep 'edge' filesep 'cell_mask'...
            filesep maskFileNames(iFrame).name];
        
        BW = imread(fileName);

        % Compute the distance transform
        D{iFrame} = single(bwdist(max(BW(:)) - BW));
        
        waitbar(iFrame / numel(maskFileNames), h);        
    end    

    % Load TM Speed map
    S = cell(numel(maskFileNames), 1);
    
    waitbar(0, h, [movieName ': Load TM speed maps...']);
    for iFile = 1:numel(speedMapFileNames{2})
        fileName = [subDirPaths{2} filesep 'post' filesep 'mat'...
            filesep speedMapFileNames{2}(iFile).name];
        load(fileName);
        [dummy1, dummy2, no] = getFilenameBody(fileName);
        no = str2double(no);
        S{no} = single(speedMap);
        clear speedMap;
        waitbar(iFile / numel(speedMapFileNames{2}), h);
    end
    % Extend the first
    fileName = [subDirPaths{2} filesep 'post' filesep 'mat'...
        filesep speedMapFileNames{2}(1).name];
    [dummy1, dummy2, no] = getFilenameBody(fileName);
    no = str2double(no);
    for i = 1:no - 1
        S{i} = S{no};
    end
    % Extend the last
    fileName = [subDirPaths{2} filesep 'post' filesep 'mat'...
        filesep speedMapFileNames{2}(end).name];
    [dummy1, dummy2, no] = getFilenameBody(fileName);
    no = str2double(no);
    for i = no+1:size(S, 1)
        S{i} = S{no};
    end
    
    % Load TM Kin map
    K = cell(numel(maskFileNames), 1);    
    
    waitbar(0, h, [movieName ': Load TM kinetics maps...']);
    for iFile = 1:numel(polyMapFileNames{2})
        fileName1 = [subDirPaths{2} filesep 'post' filesep 'mat'...
            filesep polyMapFileNames{2}(iFile).name];
        fileName2 = [subDirPaths{2} filesep 'post' filesep 'mat'...
            filesep depolyMapFileNames{2}(iFile).name];
        load(fileName1);
        load(fileName2);
        [dummy1, dummy2, no] = getFilenameBody(fileName1);
        no = str2double(no);
        K{no} = single(polyMap);
        ind = find(abs(depolyMap) > K{no});
        K{no}(ind) = single(depolyMap(ind));
        clear polyMap depolyMap;
        waitbar(iFile / numel(polyMapFileNames{2}), h);
    end
    % Extend the first
    fileName = [subDirPaths{2} filesep 'post' filesep 'mat'...
        filesep polyMapFileNames{2}(1).name];
    [dummy1, dummy2, no] = getFilenameBody(fileName);
    no = str2double(no);
    for i = 1:no - 1
        K{i} = K{no};
    end
    % Extend the last
    fileName = [subDirPaths{2} filesep 'post' filesep 'mat'...
        filesep polyMapFileNames{2}(end).name];
    [dummy1, dummy2, no] = getFilenameBody(fileName);
    no = str2double(no);
    for i = no+1:size(K, 1)
        K{i} = K{no};
    end
    
    % Gather TM information
    numTMBirths = 0;

    % tmBirths(i).loc = location of the TM birth
    % tmBirths{i}.t = birth frame
    % tmBirths(i).distToEdge = distance of TM birth to the cell edge
    % tmBirths(i).lifeTime = lifetime of TM track
    % tmBirths(i).meanSpeed = mean speed of TM track
    % tmBirths(i).meanKinScore = mean kinetics score along TM track

    array = cell(numel(find(MPM == 0)) / 2 + 1, 1);
    tmBirths = struct('loc', array, 't', array, 'distToEdge', array, ...
        'lifeTime', array, 'meanSpeed', array, 'meanKinScore', array);

    waitbar(0, h, [movieName ': Gather TM information...']);
    
    % First frame:
    waitbar(1 / numFrames, h);
    for i = 1:size(MPM, 1)
        if MPM(i, 1)
            % TM Birth
            numTMBirths = numTMBirths + 1;
            p = MPM(i, 1:2);
            tmBirths(numTMBirths).loc = p;
            tmBirths(numTMBirths).t = 1;
            tmBirths(numTMBirths).distToEdge = D{1}(p(1), p(2));

            % Gather information along the TM track
            j = 1;
            d = 0; s = 0; k = 0;
            while j < size(MPM, 2) && MPM(i, j) ~= 0
                t = ceil(j / 2);
                q = MPM(i, j:j+1);
                d = d + D{t}(q(1), q(2));
                s = s + S{t}(q(1), q(2));
                k = k + K{t}(q(1), q(2));
                j = j + 2;
            end
            tmBirths(numTMBirths).lifeTime = t;
            tmBirths(numTMBirths).meanSpeed = s / t;
            tmBirths(numTMBirths).meanKinScore = k / t;
        end
    end
    % Rest of frames:
    for j = 3:2:size(MPM, 2)
        t = ceil(j / 2);
        waitbar(t / numFrames, h);
        for i = 1:size(MPM, 1)
            if ~MPM(i, j - 1) && MPM(i, j)
                % TM Birth
                numTMBirths = numTMBirths + 1;
                p = MPM(i, j:j+1);
            tmBirths(numTMBirths).loc = p;
            tmBirths(numTMBirths).t = 1;
            tmBirths(numTMBirths).distToEdge = D{t}(p(1), p(2));

            % Gather information along the TM track
            jj = j;
            d = 0; s = 0; k = 0;
            while jj < size(MPM, 2) && MPM(i, jj) ~= 0
                tt = ceil(jj / 2);
                p = MPM(i, jj:jj+1);
                d = d + D{tt}(p(1), p(2));
                s = s + S{tt}(p(1), p(2));
                k = k + K{tt}(p(1), p(2));
                jj = jj + 2;
            end
            lifeTime = tt - t + 1;
            tmBirths(numTMBirths).lifeTime = lifeTime;
            tmBirths(numTMBirths).meanSpeed = s / lifeTime;
            tmBirths(numTMBirths).meanKinScore = k / lifeTime;
            end
        end
    end
    
    tmBirths = tmBirths(1:numTMBirths);
    
    clear MPM S K;
    
    % Load Actin MPM file
    fileName = [subDirPaths{1} filesep 'tack' filesep 'mpm.mat'];
        
    if ~exist(fileName, 'file')
        disp([movieName ': Unable to find the mpm.mat file in '...
            subDirPaths{1} '. (SKIPPING)']);
        continue;
    end
    
    load(fileName);
    
    if ~exist('MPM', 'var')
        disp([movieName ': Unable to find the Actin MPM variable in '...
            subDirPaths{1} '/tack/mpm.mat. (SKIPPING)']);
        continue;
    end
    
    clear M;
    
    % Load TM Speed map
    S = cell(numel(maskFileNames), 1);
    
    waitbar(0, h, [movieName ': Load Actin speed maps...']);
    for iFile = 1:numel(speedMapFileNames{1})
        fileName = [subDirPaths{1} filesep 'post' filesep 'mat'...
            filesep speedMapFileNames{1}(iFile).name];
        load(fileName);
        [dummy1, dummy2, no] = getFilenameBody(fileName);
        no = str2double(no);
        S{no} = single(speedMap);
        clear speedMap;
        waitbar(iFile / numel(speedMapFileNames{1}), h);
    end
    % Extend the first
    fileName = [subDirPaths{1} filesep 'post' filesep 'mat'...
        filesep speedMapFileNames{1}(1).name];
    [dummy1, dummy2, no] = getFilenameBody(fileName);
    no = str2double(no);
    for i = 1:no - 1
        S{i} = S{no};
    end
    % Extend the last
    fileName = [subDirPaths{1} filesep 'post' filesep 'mat'...
        filesep speedMapFileNames{1}(end).name];
    [dummy1, dummy2, no] = getFilenameBody(fileName);
    no = str2double(no);
    for i = no+1:size(S, 1)
        S{i} = S{no};
    end
    
    % Load Actin Kin map
    K = cell(numel(maskFileNames), 1);    
    
    waitbar(0, h, [movieName ': Load Actin kinetics maps...']);
    for iFile = 1:numel(polyMapFileNames{1})
        fileName1 = [subDirPaths{1} filesep 'post' filesep 'mat'...
            filesep polyMapFileNames{1}(iFile).name];
        fileName2 = [subDirPaths{1} filesep 'post' filesep 'mat'...
            filesep depolyMapFileNames{1}(iFile).name];
        load(fileName1);
        load(fileName2);
        [dummy1, dummy2, no] = getFilenameBody(fileName1);
        no = str2double(no);
        K{no} = single(polyMap);
        ind = find(abs(depolyMap) > K{no});
        K{no}(ind) = single(depolyMap(ind));
        clear polyMap depolyMap;
        waitbar(iFile / numel(polyMapFileNames{1}), h);
    end
    % Extend the first
    fileName = [subDirPaths{1} filesep 'post' filesep 'mat'...
        filesep polyMapFileNames{1}(1).name];
    [dummy1, dummy2, no] = getFilenameBody(fileName);
    no = str2double(no);
    for i = 1:no - 1
        K{i} = K{no};
    end
    % Extend the last
    fileName = [subDirPaths{1} filesep 'post' filesep 'mat'...
        filesep polyMapFileNames{1}(end).name];
    [dummy1, dummy2, no] = getFilenameBody(fileName);
    no = str2double(no);
    for i = no+1:size(K, 1)
        K{i} = K{no};
    end
    
    % Gather Actin information
    
    % actinEvents(i, 1) = lifetime of Actin track
    % actinEvents(i, 2) = delay from TM birth
    % actinEvents(i, 3) = mean speed of Actin track
    % actinEvents(i, 4) = mean kinetics score along Actin track
    
    % actinEvents(i, 5:6) = location of the TM birth
    % actinEvents(i, 7) = distance of TM birth to the cell edge
    % actinEvents(i, 8) = lifetime of TM track
    % actinEvents(i, 9) = mean speed of TM track
    % actinEvents(i, 10) = mean kinetics score along TM track

    numActin = 0;
    
    actinEvents = zeros(numel(find(MPM == 0)) / 2 + 1, 10);
    
    waitbar(0, h, [movieName ': Gather Actin information...']);
    for i = 1:numTMBirths
        p = tmBirths(i).loc;
        j = 2 * tmBirths(i).t - 1;
        % Find actin tracks for which distance to TM birth <= 4
        actinDist = createDistanceMatrix(p, MPM(:, j:j+1)); %#ok<COLND>
        iClosestActinTracks = find(actinDist <= 4);
        for iActinTrack = 1:numel(iClosestActinTracks)
            numActin = numActin + 1;
            ii = iClosestActinTracks(iActinTrack);
            d = 0; s = 0; k = 0;
            % Go backward on Actin track
            jmin = j;
            while jmin > 2 && MPM(ii, jmin - 2) ~= 0
                jmin = jmin - 2;
                t = ceil(jmin / 2);
                q = MPM(ii,jmin:jmin+1);
                d = d + D{t}(q(1), q(2));
                s = s + S{t}(q(1), q(2));
                k = k + K{t}(q(1), q(2));
            end
            % Go forward on Actin track
            jmax = j;
            while jmax < size(MPM, 2) && MPM(ii, jmax) ~= 0
                t = ceil(jmax / 2);
                q = MPM(ii, jmax:jmax+1);
                d = d + D{t}(q(1), q(2));
                s = s + S{t}(q(1), q(2));
                k = k + K{t}(q(1), q(2));
                jmax = jmax + 2;
            end
            lifeTime = ceil((jmax - jmin) / 2);
            
            actinEvents(numActin, 1) = lifeTime;
            actinEvents(numActin, 2) = ceil((j - jmin) / 2);
            actinEvents(numActin, 3) = s / lifeTime;
            actinEvents(numActin, 4) = k / lifeTime;
            
            actinEvents(numActin, 5:6) = tmBirths(i).loc;
            actinEvents(numActin, 7) =  tmBirths(i).distToEdge;
            actinEvents(numActin, 8) = tmBirths(i).lifeTime;
            actinEvents(numActin, 9) = tmBirths(i).meanSpeed;
            actinEvents(numActin, 10) = tmBirths(i).meanKinScore;
        end
    end
    
    actinEvents = actinEvents(1:numActin, :); %#ok<NASGU>
    
    close(h);
    
    clear MPM S K D tmBirths;

    % Save correlation scores
    save(outputFileName, 'actinEvents', 'subDirPaths');
    
    clear actinEvents;
    
    disp([movieName ': DONE']);    
end

toc;

end