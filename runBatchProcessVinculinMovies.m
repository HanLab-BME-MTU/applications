function runBatchProcessVinculinMovies

% parent directory of every movie to be analyzed
%params.rootDirectory = '/home/sb234/Projects/VinculinFA/completed/';
%params.rootDirectory = '/home/sb234/Projects/VinculinFA/completed/con/052710_con_CSUX_2';
%params.rootDirectory = '/home/sb234/Projects/VinculinFA/completed/con/062609_con_CSUX_1';
params.rootDirectory = '/Users/sylvain/Documents/Work/HMS/Projects/VinculinFA/052710_con_CSUX_2';
%params.rootDirectory = '/home/sb234/Projects/VinculinFA/SDC Y27632 20s';

% name of the channel directory subfolders
params.channelDirectory = {'ch488', 'ch560'};

% name of the processes to be run
params.procNames = {...
    'particleDetection',...
    'particleTracking',...
    'pairTracks'};
params.runSteps = [-1 -1 1];
params.batchMode = 0;

% Physical parameters
params.pixelSize = 67; %pixel size in nanometers
params.timeInterval = 10;
% N.A, magnification, camera cell size (meter), emission wavelength (meter)
sigmaPSF = getGaussianPSFsigma(1.45,100, 6.7 * 1e-6, 509 * 1e-9);

% INIT
params.setupMovieDataFunc = @setupViculinMovieData;

% PROC 1: particle detection
params.particleDetection.iChannel = 1;
params.particleDetection.detectFunc = @detectFocalAdhesionParticles;
params.particleDetection.sigmaPSF = sigmaPSF;
params.particleDetection.kSigma = 2;
params.particleDetection.alpha = .05;

% PROC 2: particle tracking
params.particleTracking.searchRadius = 5;

% PROC 3: track pairing
params.pairTracks.minOverlap = 1;
params.pairTracks.timeGap = 0;
params.pairTracks.maxEuclidianDist = 20;
params.pairTracks.sigmaPSF = sigmaPSF;
params.pairTracks.nLevels = 5;
params.pairTracks.alpha = 0.05;
params.pairTracks.probBinSize = 1e-4;

% Run all processes
batchProcessMyMovies(params);
