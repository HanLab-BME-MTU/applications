function computeLifeTime2Dist(forceRedo, varargin)
% COMPUTELIFETIME2DIST M-file for computeCorrelations.fig
%      COMPUTELIFETIME2DIST(forceRedo) compute the life time of speckle
%      track in function of its distance to the edge, formally, for each
%      Actin and TM direcotory, the readout is a matrix Nx3 where N is
%      the number of tracks and 
%      (., 1) = distance of the track birth to the edge,
%      (., 2) = distance of the track death to the edge,
%      (., 3) = life time of the track (frames),
%      'forceRedo'  forces to recompute results.

if nargin == 2
    rootDirectory = varargin{1};
else
    % Ask for the root directory.
    rootDirectory = uigetdir('', 'Select a root directory:');

    if ~ischar(rootDirectory)
        return;
    end
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

    % Check is the computation has been already performed.
    lifeTimeFileName = [path filesep 'lifeTime2Dist.mat'];    
    
    if exist(lifeTimeFileName, 'file') && ~forceRedo
        disp([movieName ': ' Param ' computation already done. (SKIPPING MOVIE)']);
        continue;
    end
        
    % Get the Actin directory name
    actinPath = [path filesep 'actin'];
    
    if ~exist(actinPath, 'dir')
        actinPath = [path filesep 'Actin'];
    end

    if ~exist(actinPath, 'dir')
        disp([movieName ': Unable to find Actin directory. (SKIPPING MOVIE)']);
        continue;
    end

    % Get the mask filename from the Actin directory.
    maskFileNames = dir([actinPath filesep 'edge' filesep...
        'cell_mask' filesep 'mask*.tif']);
    
    if isempty(maskFileNames)
        disp([movieName ': Unable to locate Actin mask files. (SKIPPING MOVIE)']);
        continue;
    end
        
    % Load the Actin MPM.
    fileName = [actinPath filesep 'tack' filesep 'mpm.mat'];
    
    if ~exist(fileName, 'file')
        disp([movieName ': Unable to find the Actin mpm.mat file. (SKIPPING MOVIE)']);
        continue;
    end
    
    load(fileName);
    
    if ~exist('MPM', 'var')
        disp([movieName ': Unable to find the Actin MPM variable. (SKIPPING MOVIE)']);
        continue;
    end
    
    MPMs = cell(2, 1);
    
    MPMs{1} = MPM;
    clear MPM;

    % Get the TM directory name
    TMPath = [path filesep 'TM2'];

    if ~exist(TMPath, 'dir')
        TMPath = [path filesep 'TM4'];
    end

    if ~exist(TMPath, 'dir')
        disp([movieName ': Unable to find TM directory. (SKIPPING MOVIE)']);
        continue;
    end

    % Load the TM MPM.
    fileName = [TMPath filesep 'tack' filesep 'mpm.mat'];
    
    if ~exist(fileName, 'file')
        disp([movieName ': Unable to find the ' TMPath ' mpm.mat file. (SKIPPING MOVIE)']);
        continue;
    end
    
    load(fileName);
    
    if ~exist('MPM', 'var')
        disp([movieName ': Unable to find the ' TMPath ' MPM variable. (SKIPPING MOVIE)']);
        continue;
    end
    
    MPMs{2} = MPM;
    clear MPM;
    
    % Get the Actin and TM tracks.
    % tracks{i} (i == 1: Actin, i == 2: TM) is a Nx6 matrix, where N is the
    % number of tracks and
    % (., 1) = Birth frame index
    % (., 2:3) = Birth position
    % (., 4) = Death frame index
    % (., 5:6) = Death position
    
    tracks = cell(2, 1);
    
    h = waitbar(0, [movieName ': Get tracks...']);
    totalTracks = size(MPMs{1}, 1) + size(MPMs{2}, 1);
    
    for i = 1:2
        for iTrack = 1:size(MPMs{i}, 1)                
            waitbar(((i - 1) * size(MPMs{1}, 1) + iTrack) / totalTracks);            
            
            j = 1;

            while j < size(MPMs{i}, 2) - 1
                birthPoint = MPMs{i}(iTrack, j:j+1);
                if birthPoint
                    birthFrame = ceil(j / 2);

                    while j < size(MPMs{i}, 2) &&...
                            MPMs{i}(iTrack, j) ~= 0
                        deathFrame = ceil(j / 2);
                        deathPoint = MPMs{i}(iTrack, j:j+1);
                        j = j + 2;                        
                    end

                    if deathFrame ~= birthFrame
                        tracks{i} = vertcat(tracks{i},...
                            [birthFrame birthPoint deathFrame deathPoint]);
                    end
                end

                j = j + 2;
            end
        end
    end
    
    close(h);
    
    % Get the frame indices from Tracks where the distance to the edge
    % needs to be computed.
    
    frameIndices = [tracks{1}(:, 1); tracks{1}(:, 4);...
        tracks{2}(:, 1); tracks{2}(:, 4)];
    frameIndices = unique(frameIndices);
    
    D = cell(numel(maskFileNames), 1);
    
    h = waitbar(0, [movieName ': Compute distance transforms...']);
    
    for i = 1:numel(frameIndices)
        waitbar(i / numel(frameIndices));
        
        iFrame = frameIndices(i);
        
        fileName = [actinPath filesep 'edge' filesep 'cell_mask'...
            filesep maskFileNames(iFrame).name];
        
        BW = imread(fileName);

        % Compute the distance transform
        D{iFrame} = bwdist(max(BW(:)) - BW);
    end
    
    close(h);
    
    % Compute lifetime
    h = waitbar(0, [movieName ': Compute lifetime...']);
    totalTracks = size(tracks{1}, 1) + size(tracks{2}, 1);
    
    R = cell(2, 1);
    R{1} = zeros(size(tracks{1}, 1), 3);
    R{2} = zeros(size(tracks{2}, 1), 3);
    
    for i = 1:2
        for iTrack = 1:size(tracks{i}, 1)
            waitbar(((i - 1) * size(tracks{1}, 1) + iTrack) / totalTracks);
            
            birthFrame = tracks{i}(iTrack, 1);
            birthPoint = tracks{i}(iTrack, 2:3);
            deathFrame = tracks{i}(iTrack, 4);
            deathPoint = tracks{i}(iTrack, 5:6);
            
            R{i}(iTrack, 1) = D{birthFrame}(birthPoint(1), birthPoint(2));
            R{i}(iTrack, 2) = D{deathFrame}(deathPoint(1), deathPoint(2));
            R{i}(iTrack, 3) = deathFrame - birthFrame + 1;
        end
    end
    
    close(h);
    
    % Save the data
    actinLifeTime2Dist = R{1}; %#ok<NASGU>
    TMLifeTime2Dist = R{2}; %#ok<NASGU>
    
    save(lifeTimeFileName, 'actinLifeTime2Dist', 'TMLifeTime2Dist');
end

toc
end