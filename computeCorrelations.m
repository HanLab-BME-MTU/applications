function computeCorrelations(param, forceRedo, varargin)
%      COMPUTECORRELATIONS(param, forceRedo) Compute the correlation
%      between Actin and TM 'param', which can be 'speed' or 'kinetics'.
%      'forceRedo' forces to recompute results.
%
%      COMPUTECORRELATIONS(param, forceRedo, rootDirectory) rootDirectory
%      is the path to movie sudirectories containing Actin and TM.

if nargin == 3
    rootDirectory = varargin{1};
else
    % Ask for the root directory.
    rootDirectory = uigetdir('', 'Select a root directory:');

    if ~ischar(rootDirectory)
        return;
    end
end

if strcmp(param, 'speed')
    paramFileName = 'speedMap*.mat';
    Param = 'Speeds';
elseif strcmp(param, 'kinetics')
    paramFileName = 'polyMap*.mat';
    Param = 'Kinetics';    
else
    error('Invalid parameters');
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

    % results format:
    %
    % corrSpeedSpeckles: [ActinSpeed TMSpeed DistanceToFront] x n
    % OR
    % corrKineticsSpeckles: [ActinKinScore TMKinScore DistanceToFront] x n
    
    corrParamFileName = [path filesep 'corr' Param '.mat'];

    % Check is the computation has been already performed.
    if exist(corrParamFileName, 'file') && ~forceRedo
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

    % Get the list of Actin cell mask file names
    maskFileNames = dir([actinPath filesep 'edge' filesep...
        'cell_mask' filesep 'mask*.tif']);
    
    if isempty(maskFileNames)
        disp([movieName ': Unable to locate Actin mask files. (SKIPPING MOVIE)']);
        continue;
    end
    
    % Get the list of Actin speckle file names.
    actinSpeckleFileNames = dir([actinPath filesep...
        'tack' filesep 'cands' filesep 'cands*.mat']);
    
    if isempty(actinSpeckleFileNames)
        disp([movieName ': Unable to locate Actin speckle files. (SKIPPING MOVIE)']);
        continue;
    end
    
    % Get the TM directory name
    TMPath = [path filesep 'TM2'];

    if ~exist(TMPath, 'dir')
        TMPath = [path filesep 'TM4'];
    end

    if ~exist(TMPath, 'dir')
        disp([movieName ': Unable to find TM directory. (SKIPPING MOVIE)']);
        continue;
    end

    % Get the list of TM speckle file names.
    TMSpeckleFileNames = dir([TMPath filesep 'tack'...
        filesep 'cands' filesep 'cands*.mat']);
    
    if isempty(TMSpeckleFileNames)
        disp([movieName ': Unable to locate TM speckle files. (SKIPPING MOVIE)']);
        continue;
    end

    % Get the list of Actin speed|kinetics map file names.
    actinParamFileNames = dir([actinPath filesep...
        'post' filesep 'mat' filesep paramFileName]);
            
    if isempty(actinParamFileNames)
        disp([movieName ': Unable to locate Actin ' param ' map files. (SKIPPING MOVIE)']);
        continue;
    end
            
    % Get the list of TM speed|kinetics map file names.
    TMParamFileNames = dir([TMPath filesep 'post'...
        filesep 'mat' filesep paramFileName]);
    
    % Check the number of param map in both Actin and TM are equal.
    if numel(actinParamFileNames) ~= numel(TMParamFileNames)
        disp([movieName ': Number of Actin and TM ' param ' map files differ. (SKIPPING MOVIE)']);
        continue;
    end
                
    numFrames = numel(actinParamFileNames);

    corrParam = [];
    
    h = waitbar(0, movieName);
    
    for iFile = 1:numFrames
        frameName = ['Frame ' num2str(iFile) '/' num2str(numFrames)];
        
        waitbar(iFile / numFrames, h, [movieName ' - ' frameName...
            ': processing...']);
        
        % Get the index of param map file.
        [dummy1, dummy2, no] = getFilenameBody(...
            actinParamFileNames(iFile).name);
        no = str2double(no);
        
        % Load the Actin param map.
        fileName = [actinPath filesep 'post' filesep 'mat'...
            filesep actinParamFileNames(iFile).name];
        load(fileName);
        if strcmp(param, 'speed')            
            ActinParamMap = speedMap;
            clear speedMap;
        else
            ActinParamMap = polyMap;
            [path2, body2, no2, ext2] = getFilenameBody(fileName);
            load([path2 filesep 'depolyMap_' no2 ext2]);
            ind = find(abs(depolyMap) > ActinParamMap);
            ActinParamMap(ind) = depolyMap(ind);
            clear polyMap depolyMap;
        end
        
        % Load the TM speed map.
        fileName = [TMPath filesep 'post' filesep 'mat' filesep...
            TMParamFileNames(iFile).name];
        load(fileName);
        if strcmp(param, 'speed')
            TMParamMap = speedMap;
            clear speedMap;
        else
            TMParamMap = polyMap;
            [path2, body2, no2, ext2] = getFilenameBody(fileName);
            load([path2 filesep 'depolyMap_' no2 ext2]);
            ind = find(abs(depolyMap) > TMParamMap);
            TMParamMap(ind) = depolyMap(ind);
            clear polyMap depolyMap;
        end
        
        % Load the mask.
        fileName = [actinPath filesep 'edge' filesep 'cell_mask'...
            filesep maskFileNames(no).name];        
        BW = imread(fileName);
        
        % Load the Actin speckle.
        fileName = [actinPath filesep 'tack' filesep 'cands'...
            filesep actinSpeckleFileNames(no).name];        
        load(fileName);
        status = cat(1, cands(:).status); %#ok<COLND>
        ActinSpeckles = cat(1, cands(:).Lmax); %#ok<COLND>
        ActinSpeckles = ActinSpeckles(status == 1, :);
        clear cands;

        % Load the TM speckle.
        fileName = [TMPath filesep 'tack' filesep 'cands'...
            filesep actinSpeckleFileNames(no).name];
        load(fileName);
        status = cat(1, cands(:).status); %#ok<COLND>
        TMSpeckles = cat(1, cands(:).Lmax); %#ok<COLND>
        TMSpeckles = TMSpeckles(status == 1, :);
        clear cands;

        % Compute the distance transform
        D = bwdist(max(BW(:)) - BW);

        ind1 = sub2ind(size(D), ActinSpeckles(:, 1), ActinSpeckles(:, 2));
        ind2 = sub2ind(size(D), TMSpeckles(:, 1), TMSpeckles(:, 2));
        
        % Find speckles indices belonging to the narrow band
        %ind11 = ind1(loDist <= D(ind1) & D(ind1) <= hiDist);
        %ind22 = ind2(loDist <= D(ind2) & D(ind2) <= hiDist);
                
        % Get the Actin and TM param value at every speckle position.
        corrParam = vertcat(corrParam,...
            horzcat(ActinParamMap(ind1), TMParamMap(ind1), D(ind1)),...
            horzcat(ActinParamMap(ind2), TMParamMap(ind2), D(ind2)));
    end
    
    % Save total corrSpeeds
    save(corrParamFileName, 'corrParam');
    
    close(h);
    
    disp([movieName ': DONE']);    
end

toc

end