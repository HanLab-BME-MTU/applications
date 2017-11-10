% List all TIFF files in the main folder
% if strcmp(getenv('USER'),'kj35')
mainFolder = 'C:\kjData\maki\newUnarchived\JulieSebastien\1306_cambridgeData\metaMovies';
% elseif strcmp(getenv('USER'),'sebastien')
%     mainFolder = fullfile(getenv('HOME'),'Documents','Julie','testMovie15');
% end

tiffFiles = dir(fullfile(mainFolder,'*.tif*'));
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
        MD(i).getPackage(1).createDefaultProcess(1);
        MD(i).getPackage(1).createDefaultProcess(2);
        
        % Create spindle axis and sister grouping processes
        MD(i).addProcess(SpindlePolesEBProcess(MD(i)));
        MD(i).addProcess(SisterGroupingProcess(MD(i)));
        
        % Create k-EB detection processes
        MD(i).addProcess(KEBDetectionProcess(MD(i)));
        
    end
    
end

%% Set analysis parameters

%NOTE: space units = pixels
%       time units = frames

for i=1:nMovies
    
    % Gaussian mixture-model fitting
    %general
    funParams = MD(i).getProcess(1).funParams_;
    funParams.ChannelIndex=2; % Detect mCherry-CENPA objects
    %function-specific
    funParams.verbose = 0;
    funParams.detectionParam.psfSigma = 1.9;
    funParams.detectionParam.bitDepth = 16;
    funParams.detectionParam.alphaLocMax = 0.1;
    funParams.detectionParam.integWindow = 0;
    funParams.detectionParam.doMMF = 1;
    funParams.detectionParam.testAlpha = struct('alphaR',0.0001,'alphaA',0.1,'alphaD',0.05,'alphaF',0);
    funParams.detectionParam.numSigmaIter = 0;
    funParams.detectionParam.visual = 0;
    funParams.detectionParam.background = [];
    %general
    parseProcessParams(MD(i).getProcess(1),funParams);
    
    % Tracking
    %general
    funParams = MD(i).getProcess(2).funParams_;
    funParams.ChannelIndex=2; % Track mCherry-CENPA objects
    %function-specific
    funParams.verbose = 0;
    %gap closing
    funParams.gapCloseParam.timeWindow = 6;
    funParams.gapCloseParam.mergeSplit = 0;
    funParams.gapCloseParam.minTrackLen = 1;
    funParams.gapCloseParam.diagnostics = [];
    %cost matrix 1
    funParams.costMatrices(1).parameters.linearMotion = 0;
    funParams.costMatrices(1).parameters.minSearchRadius = 2;
    funParams.costMatrices(1).parameters.maxSearchRadius = 5;
    funParams.costMatrices(1).parameters.brownStdMult = 3;
    funParams.costMatrices(1).parameters.useLocalDensity = 1;
    funParams.costMatrices(1).parameters.nnWindow = funParams.gapCloseParam.timeWindow;
    funParams.costMatrices(1).parameters.kalmanInitParam = [];
    funParams.costMatrices(1).parameters.diagnostics = [];
    %cost matrix 2
    funParams.costMatrices(2).parameters.linearMotion = 0;
    funParams.costMatrices(2).parameters.minSearchRadius = 2;
    funParams.costMatrices(2).parameters.maxSearchRadius = 5;
    funParams.costMatrices(2).parameters.brownStdMult = 3*ones(funParams.gapCloseParam.timeWindow,1);
    funParams.costMatrices(2).parameters.brownScaling = [0.25 0.01];
    funParams.costMatrices(2).parameters.timeReachConfB = funParams.gapCloseParam.timeWindow;
    funParams.costMatrices(2).parameters.ampRatioLimit = [0.7 4];
    funParams.costMatrices(2).parameters.lenForClassify = 5;
    funParams.costMatrices(2).parameters.useLocalDensity = 0;
    funParams.costMatrices(2).parameters.nnWindow = funParams.gapCloseParam.timeWindow;
    funParams.costMatrices(2).parameters.linStdMult = 1*ones(funParams.gapCloseParam.timeWindow,1);
    funParams.costMatrices(2).parameters.linScaling = [0.25 0.01];
    funParams.costMatrices(2).parameters.timeReachConfL = funParams.gapCloseParam.timeWindow;
    funParams.costMatrices(2).parameters.maxAngleVV = 30;
    funParams.costMatrices(2).parameters.gapPenalty = 1.5;
    funParams.costMatrices(2).parameters.resLimit = [];
    %general
    parseProcessParams(MD(i).getProcess(2),funParams);
    
    % Spindle poles
    %general
    funParams = MD(i).getProcess(3).funParams_;
    funParams.ChannelIndex=1; % Derive spindle axis from GFP-EB3 images
    %function-specific
    funParams.doPlot = 0;  
    funParams.numPoles = 2;
    %general
    parseProcessParams(MD(i).getProcess(3), funParams);
    
    % Sister pairing
    %general
    funParams = MD(i).getProcess(4).funParams_;
    funParams.ChannelIndex=2; % Group mCherry-CENPA tracks
    %function-specific
    funParams.maxAngle = 45*pi/180; %radians
    funParams.maxDist = 20; %pixels
    funParams.minOverlap = 10; %frames
    funParams.useAlignment = 1;
    funParams.robust = 0;
    funParams.associateSis2Pole = 1;
    %general
    parseProcessParams(MD(i).getProcess(4), funParams);
    
    % k-EB detection
    %general
    funParams = MD(i).getProcess(5).funParams_;
    funParams.ChannelIndex=1; % Detect GFP-EB3 signal
    %function-specific
    funParams.radiusEB = 3; %pixels
    funParams.lengthAlongMT = 7;
    %general
    parseProcessParams(MD(i).getProcess(5), funParams);

end

%% Run analysis processes
for i=1:nMovies
    cellfun(@run,MD(i).processes_);
end

%% Set visualization options
for i=1:nMovies
    % Get default display class for spindle axis/spindle poles
    poleDisplayMethod = MD(i).getProcess(3).getDrawableOutput.defaultDisplayMethod(1);
    poleDisplayMethod.Color = [1 1 0]; % Set color to yellow
    MD(i).getProcess(3).setDisplayMethod(1, 1, poleDisplayMethod);
    % Get default display class for sister grouping
    pairDisplayMethod = MD(i).getProcess(4).getDrawableOutput.defaultDisplayMethod(2);
    pairDisplayMethod.Color = [1 1 0]; % Set color to yellow
    MD(i).getProcess(4).setDisplayMethod(1, 2, pairDisplayMethod);
end
