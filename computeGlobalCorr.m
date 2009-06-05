function computeGlobalCorr(varargin)
%      COMPUTEGLOBALCORR(forceRedo) Compute the correlation
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
    outputFileName = [path filesep 'globalCorrelations.mat'];
    
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
%     polyMapFileNames = cell(2, 1);
%     depolyMapFileNames = cell(2, 1);
    
    error = 0;
    
    for k = 1:2
        speedMapFileNames{k} = dir([subDirPaths{k} filesep...
            'post' filesep 'mat' filesep 'speedMap*.mat']);
        
        if isempty(speedMapFileNames{k})
            disp([movieName ': Unable to locate speed map files. (SKIPPING)']);
            error = 1;
            break;
        end
        
%         polyMapFileNames{k} = dir([subDirPaths{k} filesep...
%             'post' filesep 'mat' filesep 'polyMap*.mat']);
%         
%         if isempty(polyMapFileNames{k})
%             disp([movieName ': Unable to locate poly map files. (SKIPPING)']);
%             error = 1;
%             break;
%         end
%         
%         depolyMapFileNames{k} = dir([subDirPaths{k} filesep...
%             'post' filesep 'mat' filesep 'depolyMap*.mat']);
%         
%         if isempty(polyMapFileNames{k})
%             disp([movieName ': Unable to locate depoly map files. (SKIPPING)']);
%             error = 1;
%             break;
%         end
    end
    
    if error
        continue;
    end

    % Load MPM file
    MPMs = cell(2, 1);
    
    for k = 1:2
        fileName = [subDirPaths{k} filesep 'tack' filesep 'mpm.mat'];
        
        if ~exist(fileName, 'file')
            disp([movieName ': Unable to find the mpm.mat file in '...
                subDirPaths{k} '. (SKIPPING)']);
            error = 1;
            break;
        end
    
        load(fileName);
    
        if ~exist('MPM', 'var')
            disp([movieName ': Unable to find the Actin MPM variable in '...
                subDirPaths{k} '/tack/mpm.mat. (SKIPPING)']);
            error = 1;
            break;
        end

        MPMs{k} = MPM;
        
        clear M MPM;
    end

    if error
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
        D{iFrame} = single(bwdist(max(BW(:)) - BW));
        
        waitbar(iFrame / numel(maskFileNames), h);        
    end    
    
    % Compute track mean distance to the edge
    %
    % Format:
    % R(i, 1:2) = mean position of the track i
    % R(i, 3) = mean distance of the track i to the cell edge
    % R(i, 4) = mean actin speed along track i
    % R(i, 5) = mean TM speed score along track i
    % R(i, 6) = life time of track i  
    
    R = cell(2, 1);
    for k = 1:2
        numTracks = 0;
        R{k} = zeros(numel(find(MPMs{k} == 0)) / 2 + 1, 6, 'single');
        
        waitbar(0, h, [movieName ': Compute track mean distance to the edge...']);
        for i = 1:size(MPMs{k}, 1)                
            waitbar(i / size(MPMs{k}, 1), h);            
            
            j = 1;

            while j < size(MPMs{k}, 2) - 1
                if MPMs{k}(i, j:j+1)
                    birth = ceil(j / 2);
                    numTracks = numTracks + 1;
                    d = 0;
                    m = zeros(1, 2);
                    while j < size(MPMs{k}, 2) && MPMs{k}(i, j) ~= 0
                        iFrame = ceil(j / 2);
                        p = MPMs{k}(i, j:j+1);
                        death = iFrame;
                        d = d + D{iFrame}(p(1), p(2));
                        m = m + p;
                        j = j + 2;                        
                    end
                    lifeTime = death - birth + 1;
                    R{k}(numTracks, 1:2) = m / lifeTime;
                    R{k}(numTracks, 3) = d / lifeTime;
                    R{k}(numTracks, 6) = lifeTime;
                end

                j = j + 2;
            end
        end
        
        R{k} = R{k}(1:numTracks, :);        
    end    
    
    clear D;
    
    % Load the speed maps
    S = cell(numel(maskFileNames), 2);
    
    for k = 1:2
        waitbar(0, h, [movieName ': Load speed maps (' num2str(k) '/2)...']);
        for iFile = 1:numel(speedMapFileNames{k})
            fileName = [subDirPaths{k} filesep 'post' filesep 'mat'...
                filesep speedMapFileNames{k}(iFile).name];
            load(fileName);
            [dummy1, dummy2, no] = getFilenameBody(fileName);
            no = str2double(no);
            S{no, k} = single(speedMap);
            clear speedMap;
            waitbar(iFile / numel(speedMapFileNames{k}), h);
        end
        % Extend the first
        fileName = [subDirPaths{k} filesep 'post' filesep 'mat'...
            filesep speedMapFileNames{k}(1).name];
        [dummy1, dummy2, no] = getFilenameBody(fileName);
        no = str2double(no);
        for i = 1:no - 1
            S{i, k} = S{no, k};
        end
        % Extend the last
        fileName = [subDirPaths{k} filesep 'post' filesep 'mat'...
            filesep speedMapFileNames{k}(end).name];
        [dummy1, dummy2, no] = getFilenameBody(fileName);
        no = str2double(no);
        for i = no+1:size(S, 1)
            S{i, k} = S{no, k};
        end
    end
    
    % Compute track mean speed
    for k = 1:2
        numTracks = 0;

        waitbar(0, h, [movieName ': Compute track mean speed (' num2str(k) '/2)...']);
        for i = 1:size(MPMs{k}, 1)                
            waitbar(i / size(MPMs{k}, 1), h);            
            
            j = 1;

            while j < size(MPMs{k}, 2) - 1
                if MPMs{k}(i, j:j+1)
                    birth = ceil(j / 2);
                    numTracks = numTracks + 1;
                    s1 = 0; s2 = 0;
                    while j < size(MPMs{k}, 2) && MPMs{k}(i, j) ~= 0
                        iFrame = ceil(j / 2);
                        p = MPMs{k}(i, j:j+1);
                        death = iFrame;
                        s1 = s1 + S{iFrame, 1}(p(1), p(2));
                        s2 = s2 + S{iFrame, 2}(p(1), p(2));
                        j = j + 2;                        
                    end
                    lifeTime = death - birth + 1;
                    R{k}(numTracks, 4) = s1 / lifeTime;
                    R{k}(numTracks, 5) = s2 / lifeTime;
                end

                j = j + 2;
            end
        end
    end
    
    clear S;

%     % Load the kinetics maps
%     K = cell(numel(maskFileNames), 2);    
%     
%     for k = 1:2        
%         waitbar(0, h, [movieName ': Load kinetics maps (' num2str(k) '/2)...']);
%         for iFile = 1:numel(polyMapFileNames{k})
%             fileName1 = [subDirPaths{k} filesep 'post' filesep 'mat'...
%                 filesep polyMapFileNames{k}(iFile).name];
%             fileName2 = [subDirPaths{k} filesep 'post' filesep 'mat'...
%                 filesep depolyMapFileNames{k}(iFile).name];
%             load(fileName1);
%             load(fileName2);
%             [dummy1, dummy2, no] = getFilenameBody(fileName1);
%             no = str2double(no);
%             K{no, k} = single(polyMap);
%             ind = find(abs(depolyMap) > K{no, k});
%             K{no, k}(ind) = single(depolyMap(ind));
%             clear polyMap depolyMap;
%             waitbar(iFile / numel(speedMapFileNames{k}), h);
%         end
%         % Extend the first
%         fileName = [subDirPaths{k} filesep 'post' filesep 'mat'...
%             filesep polyMapFileNames{k}(1).name];
%         [dummy1, dummy2, no] = getFilenameBody(fileName);
%         no = str2double(no);
%         for i = 1:no - 1
%             K{i, k} = K{no, k};
%         end
%         % Extend the last
%         fileName = [subDirPaths{k} filesep 'post' filesep 'mat'...
%             filesep polyMapFileNames{k}(end).name];
%         [dummy1, dummy2, no] = getFilenameBody(fileName);
%         no = str2double(no);
%         for i = no+1:size(K, 1)
%             K{i, k} = K{no, k};
%         end        
%     end
% 
%     % Compute track mean kin score
%     for k = 1:2
%         numTracks = 0;
%         
%         waitbar(0, h, [movieName ': Compute track mean kin score (' num2str(k) '/2)...']);
%         for i = 1:size(MPMs{k}, 1)                
%             waitbar(i / size(MPMs{k}, 1), h);            
%             
%             j = 1;
% 
%             while j < size(MPMs{k}, 2) - 1
%                 if MPMs{k}(i, j:j+1)
%                     birth = ceil(j / 2);
%                     numTracks = numTracks + 1;
%                     k1 = 0; k2 = 0;
%                     while j < size(MPMs{k}, 2) && MPMs{k}(i, j) ~= 0
%                         iFrame = ceil(j / 2);
%                         p = MPMs{k}(i, j:j+1);
%                         death = iFrame;
%                         k1 = k1 + K{iFrame, 1}(p(1), p(2));
%                         k2 = k2 + K{iFrame, 2}(p(1), p(2));                        
%                         j = j + 2;                        
%                     end
%                     lifeTime = death - birth + 1;
%                     R{k}(numTracks, 4) = k1 / lifeTime;
%                     R{k}(numTracks, 5) = k2 / lifeTime;
%                 end
% 
%                 j = j + 2;
%             end
%         end
%     end
% 
%     clear K;
    
    % Save correlation scores
    save(outputFileName, 'R', 'subDirPaths');
    
    clear R MPMs;
    
    close(h);
    
    disp([movieName ': DONE']);    
end

toc;

end