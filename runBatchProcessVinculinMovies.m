function runBatchProcessVinculinMovies

% parent directory of every movie to be analyzed
params.rootDirectory = '/home/sb234/Projects/VinculinFA/completed/con';
%params.rootDirectory = '/home/sb234/Projects/VinculinFA/completed/cre/102709_cre_CSUX_8';
%params.rootDirectory = '/home/sb234/Projects/VinculinFA/completed/con/052710_con_CSUX_2';
%params.rootDirectory = '/home/sb234/Projects/VinculinFA/completed/con/052710_con_CSUX_5';
%params.rootDirectory = '/home/sb234/Projects/VinculinFA/SDC Y27632 20s';

% name of the channel directory subfolders
params.channelDirectory = {'ch488', 'ch560'};

% name of the processes to be run
params.procNames = {...
    'distanceTransform',...
    'particleDetection',...
    'particleTracking',...
    'pairTracks',...
    'contours',...
    'protrusion',...
    'windows',...
    'protrusionSamples',...
    'labels',...
    'figures'};

params.runSteps = [0 0 0 0 0 0 0 0 0 1];
params.batchMode = false;

% Physical parameters
params.pixelSize = 67; %pixel size in nanometers
params.timeInterval = 10;
% N.A, magnification, camera cell size (meter), emission wavelength (meter)
sigmaPSF = getGaussianPSFsigma(1.45,100, 6.7 * 1e-6, 509 * 1e-9);

% INIT
params.setupMovieDataFunc = @setupViculinMovieData;

% PROC 1: distance transform
params.distanceTransform.required = struct();
params.distanceTransform.optional = struct();

% PROC 2: particle detection
params.particleDetection.required.iChannel = 1;
params.particleDetection.required.sigmaPSF = sigmaPSF;
params.particleDetection.optional.mode = 'xyArtc';
params.particleDetection.optional.kSigma = 3;
params.particleDetection.optional.alpha = .01;
params.particleDetection.optional.minDist = .25;

% PROC 3: particle tracking
params.particleTracking.required.searchRadius = 5;
params.particleTracking.optional = struct();

% PROC 4: track pairing
params.pairTracks.required = struct();
params.pairTracks.optional.minLifetime = 3;
params.pairTracks.optional.maxDistance = 1675; % nm (25 pixels)
params.pairTracks.optional.minOverlap = 1;
params.pairTracks.optional.bandWidth = 1000;   % nm
params.pairTracks.optional.minDistance = 335;  % nm (5 pixels)
params.pairTracks.optional.alpha = 0.05;

% PROC 5: contours
dContour = 1000 / params.pixelSize; % 1 um
params.contours.required.distVals = 0:dContour:500;
params.contours.optional.forceClose = false;
params.contours.optional.maskChannels = 1;
params.contours.optional.contourName = ['contours_'  num2str(dContour) 'pix.mat'];
params.contours.optional.contourAlignAlg = 1;

% PROC 6: protrusion
params.protrusion.required = struct();
params.protrusion.optional.maskChannel = 1;
params.protrusion.optional.downSample = 20;
params.protrusion.optional.nSeg = 30;

% PROC 7: windows
params.windows.required = struct();
params.windows.optional.methodStr = 'p';
params.windows.optional.winSize = 1000 / params.pixelSize; % ~1um;
params.windows.optional.iOuter = 2;
params.windows.optional.iInner = 4;
params.windows.optional.windowName = ['windows_' ...
    params.windows.optional.methodStr '_' ...
    num2str(dContour) 'by' num2str(params.windows.optional.winSize) ...
    'pix_' num2str(params.windows.optional.iOuter) '_' ...
    num2str(params.windows.optional.iInner) '.mat'];

% PROC 8: protrusion sampling
params.protrusionSamples.required = struct();
params.protrusionSamples.optional.protName = ['protSamples_' ...
    params.windows.optional.methodStr '_' ...
    params.windows.optional.windowName  '.mat'];

% PROC 9: window labeling
params.labels.required.method = 'sector';
params.labels.optional = struct();

% PROC 10: figures
params.figures.required.bandWidth = 12000;     % nm
params.figures.optional.minActinLifetime = 3;
params.figures.optional.minSegmentsPerBin = 15;
params.figures.optional.alpha = 1e-5;

% Run all processes
batchProcessMyMovies(params);
