function runBatchProcessVinculinMovies

% parent directory of every movie to be analyzed
%params.rootDirectory = '/home/sb234/Projects/VinculinFA/completed/con/052710_con_CSUX_2';
params.rootDirectory = '/home/sb234/Projects/VinculinFA/completed/cre/060910_cre_CSUX_1';
%params.rootDirectory = '/home/sb234/Projects/VinculinFA/completed/cre/060910_cre_CSUX_1';
%params.rootDirectory = '/Users/sylvain/Documents/Work/HMS/Projects/VinculinFA/062309_cre_CSUX_3';
%params.rootDirectory = '/home/sb234/Projects/VinculinFA/SDC Y27632 20s';

% name of the channel directory subfolders
params.channelDirectory = {'ch488', 'ch560'};

% name of the processes to be run
params.procNames = {...
    'distanceTransform',...
    'particleDetection',...
    'particleTracking',...
    'pairTracks',...
    'figures'};

params.runSteps = [0 0 0 1 1];
params.batchMode = 0;

% Physical parameters
params.pixelSize = 67; %pixel size in nanometers
params.timeInterval = 10;
% N.A, magnification, camera cell size (meter), emission wavelength (meter)
sigmaPSF = getGaussianPSFsigma(1.45,100, 6.7 * 1e-6, 509 * 1e-9);

% INIT
params.setupMovieDataFunc = @setupViculinMovieData;

% PROC 1: distance transform
params.distanceTransform = struct();

% PROC 2: particle detection
params.particleDetection.iChannel = 1;
params.particleDetection.detectFunc = @detectFocalAdhesionParticles;
params.particleDetection.sigmaPSF = sigmaPSF;
params.particleDetection.kSigma = 3;
params.particleDetection.alpha = .05;
params.particleDetection.minDist = .2;

% PROC 3: particle tracking
params.particleTracking.searchRadius = 5;

% PROC 4: track pairing

params.pairTracks.maxDistance = 2010; % nm (30 pixels)
params.pairTracks.minOverlap = 1;
params.pairTracks.bandWidth = 1000;   % nm
params.pairTracks.minDistance = 335;  % nm (5 pixels)
params.pairTracks.alpha = 0.05;

% PROC 5: figures
params.figures = struct();

% Run all processes
batchProcessMyMovies(params);
