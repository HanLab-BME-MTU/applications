% List all TIFF files in the main folder
mainFolder = '/Users/sebastien/Documents/Julie/kMTProject';
tiffFiles = dir(fullfile(mainFolder,'*.tif'));
nMovies = numel(tiffFiles);

% Create array of movies
MD(nMovies,1) = MovieData();
for i=1:nMovies
    % Create MovieData object using Bioformats
    % For a moviexxxx.tif file, a moviexx folder will be created where the
    % moviexxx.mat file and all analysis will be stored
    MD(i) = MovieData.load(fullfile(mainFolder,tiffFiles(i).name));
end

% Create list of movies
listFolder = fullfile(mainFolder,'kMTAnalysis');
if ~isdir(listFolder),  mkdir(listFolder); end
ML=MovieList(MD,listFolder);
ML.setPath(listFolder);
ML.setFilename('kMTList.mat');
ML.save;


%% Recreate analysis processes

% Flag to reset all analysis
resetAnalysis = true; 

if resetAnalysis
    arrayfun(@reset,MD); % Wipe out all existing analysis
    for i=1:nMovies
        % Create UTrackPackage and its first 2 processes (detection &
        % tracking)
        MD(i).addPackage(UTrackPackage(MD(i)));
        MD(i).packages_{1}.createDefaultProcess(1);
        MD(i).packages_{1}.createDefaultProcess(2);
        
        % Create sister grouping and kMT detection processes
        MD(i).addProcess(SisterGroupingProcess(MD(i)));
        MD(i).addProcess(KMTDetectionProcess(MD(i)));
    end
end

%% Set analysis parameters
for i=1:nMovies
    
    % Set parameters for the Gaussian mixture-model fitting
    funParams = MD(i).processes_{1}.funParams_;
    funParams.ChannelIndex=2; % Detect mCherry-CENPA objects
    funParams.detectionParam.psfSigma=1.5; % Set the psfSigma
    parseProcessParams(MD(i).processes_{1},funParams);
    
    % Set tracking parameters
    funParams = MD(i).processes_{2}.funParams_;
    funParams.ChannelIndex=2; % Track mCherry-CENPA objects
    funParams.costMatrices(1).parameters.diagnostics=[]; % Do not output diagnostics
    parseProcessParams(MD(i).processes_{2},funParams);
    
    % Set sister pairing parameters
    funParams = MD(i).processes_{3}.funParams_;
    funParams.ChannelIndex=2; % Group mCherry-CENPA tracks 
    parseProcessParams(MD(i).processes_{3}, funParams);
    
    % Set kMT detection parameters
    funParams = MD(i).processes_{4}.funParams_;
    funParams.ChannelIndex=1; % Detect gfp-EB3 comets
    parseProcessParams(MD(i).processes_{4}, funParams);

end

%% Run analysis processes
for i=1:nMovies
    cellfun(@run,MD(i).processes_);
end
