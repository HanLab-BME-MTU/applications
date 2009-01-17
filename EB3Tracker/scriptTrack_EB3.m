
%% general gap closing parameters

% maximum allowed time gap (in frames) between a track segment end and a
% track segment start that allows linking them. if timeWindow = n, then
% there are n-1 frames between the last detection event of the first track
% segment and the first detection event of the second track segment.
gapCloseParam.timeWindow = 8; 

% 1 if merging and splitting are to be considered, 0 otherwise.
% even though EB comets might overlap in a given frame, they do not merge
% or split physically, so this should be zero.
gapCloseParam.mergeSplit = 0;

% minimum length of track segments (in frames) from linking to be used in
% gap closing. 1 if everything, 2 if all tracks 2 or more frames long, etc.
gapCloseParam.minTrackLen = 3; 

%% cost matrix for frame-to-frame linking

%function name
costMatrices(1).funcName = 'costMatLinearMotionLink_EB3'; 

%parameters

parameters.linearMotion = 1; %use linear motion Kalman filter.

parameters.minSearchRadius = 10; %minimum allowed search radius. The search radius is calculated on the spot in the code given a feature's motion parameters. If it happens to be smaller than this minimum, it will be increased to the minimum.
parameters.maxSearchRadius = 15; %maximum allowed search radius. Again, if a feature's calculated search radius is larger than this maximum, it will be reduced to this maximum.
parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.

parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).


parameters.kalmanInitParam.searchRadiusFirstIteration = 20; %Kalman filter initialization parameters.
%parameters.kalmanInitParam = [];

costMatrices(1).parameters = parameters;
clear parameters

%% cost matrix for gap closing

%function name
costMatrices(2).funcName = 'costMatLinearMotionCloseGaps_EB3';

%parameters

%needed all the time
% parameters.minSearchRadius = 1; %minimum allowed search radius.
% parameters.maxSearchRadius = 10; %maximum allowed search radius.
% 
% parameters.brownStdMult = 3*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate Brownian search radius from standard deviation.
% parameters.timeReachConfB = 2; %in the code, the search radius expands with the time gap (since a particle is expected to move further away in a longer gap than in a shorter one). This parameter controls how fast the search radius grows with time. timeReachConfB stands for time to reach confinement for the Brownian part of the motion. So before timeReachConfB, the search radius grows with the square root of time, after that it grows very, very slowly (it's almost fixed).
% parameters.linStdMult = 1*ones(gapCloseParam.timeWindow,1); %multiplication factor to calculate linear search radius from standard deviation.
% parameters.timeReachConfL = gapCloseParam.timeWindow; %same as timeReachConfB, but for the linear part of the motion.
% 
% parameters.lenForClassify = 5; %minimum track segment length to classify it as linear or random.
% 
% parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the gap closing and merging/splitting step.
parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before/after the current one where you want to look for a track's nearest neighbor at its end/start (in the gap closing step).

parameters.maxFAngle = 45;
parameters.maxBAngle = 15;
parameters.backVelMultFactor = 2;

costMatrices(2).parameters = parameters;
clear parameters

%% Kalman filter function names

kalmanFunctions.reserveMem = 'kalmanResMemLM';
kalmanFunctions.initialize = 'kalmanInitLinearMotion_EB3';
kalmanFunctions.calcGain = 'kalmanGainLinearMotion_EB3';
kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';

%% additional input

%saveResults
trackDir = [pwd filesep 'track'];
if ~isdir(trackDir)
    mkdir(trackDir)
end
saveResults.dir = trackDir;  %directory where to save input and output
saveResults.filename = 'trackResults.mat'; %name of file where input and output are saved


%verbose
verbose = 1;

%problem dimension
probDim = 2;

%% tracking function call

[tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalman(movieInfo,...
    costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

