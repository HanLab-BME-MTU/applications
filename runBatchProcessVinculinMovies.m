function runBatchProcessVinculinMovies

% parent directory of every movie to be analyzed
%params.rootDirectory = '/home/sb234/Projects/VinculinFA/completed/';
params.rootDirectory = '/Users/sylvain/Documents/Work/HMS/Projects/VinculinFA/052710_con_CSUX_2';

% name of the channel directory subfolders
params.channelDirectory = {'ch488', 'ch560'};

% name of the processes to be run
params.procNames = {...
    'particleDetection',...
    'particleTracking',...
    'pairTracks'};
params.runSteps = [0 1 -1];
params.batchMode = 1;

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

% PROC 2: particle tracking
params.particleTracking.searchRadius = 5;

% PROC 3: vimentin track pairing
params.pairTracks.iChannel = 1;
params.pairTracks.sigmaPSF = sigmaPSF;
params.pairTracks.thetaTh = pi/16;
params.pairTracks.alpha = 0.05;

% Run all processes
batchProcessMyMovies(params);
