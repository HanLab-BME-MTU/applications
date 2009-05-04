function plusTipTracker(runInfo,timeWindow,minTrackLen,minRadius,maxRadius,maxFAngle,maxBDist)
% plusTipTracker is the tracking function

%% get runInfo in correct format
if nargin<1 || isempty(runInfo)
    homeDir=pwd;
    runInfo.anDir=uigetdir(pwd,'Please select analysis directory');
    cd([runInfo.anDir filesep '..'])
    runInfo.imDir=uigetdir(pwd,'Please select image directory');
    cd(homeDir)
else
    % adjust for OS
    if ~isfield(runInfo,'imDir') || ~isfield(runInfo,'anDir')
        error('--plusTipTracker: first argument should be a structure with fields imDir and anDir');
    else
        [runInfo.anDir] = formatPath(runInfo.anDir);
        homeDir=pwd;
        cd(runInfo.anDir)
        [runInfo.imDir] = formatPath(runInfo.imDir);
        cd(homeDir)
    end
end

% load movieInfo (detection result)
featDir  = [runInfo.anDir filesep 'feat'];
if ~isdir(featDir)
    error('--plusTipTracker: feat directory missing')
else
    if exist([featDir filesep 'movieInfo.mat'],'file')
        load([featDir filesep 'movieInfo.mat'])
    else
        error('--plusTipTracker: movieInfo missing...')
    end
end

%% general gap closing parameters
% maximum allowed time gap (in frames) between a track segment end and a
% track segment start that allows linking them. if timeWindow = n, then
% there are n-1 frames between the last detection event of the first track
% segment and the first detection event of the second track segment.
if nargin<2 || isempty(timeWindow)
    error('--plusTipTracker: time window input missing')
else
    gapCloseParam.timeWindow = timeWindow;
end
% minimum length of track segments (in frames) from linking to be used in
% gap closing. 1 if everything, 2 if all tracks 2 or more frames long, etc.
if nargin<3 || isempty(minTrackLen)
    error('--plusTipTracker: min track length input missing')
else
    gapCloseParam.minTrackLen = minTrackLen;
end

% 1 if merging and splitting are to be considered, 0 otherwise.
% even though EB comets might overlap in a given frame, they do not merge
% or split physically, so this should be zero.
gapCloseParam.mergeSplit = 0;

%% cost matrix for frame-to-frame linking

%function name
costMatrices(1).funcName = 'costMatLinearMotionLink_EB3'; 

%used 10 and 15 for 2sec data,3/5 for Claudio 
%minimum allowed search radius. the search radius is calculated on the spot
%in the code given a feature's motion parameters. if it happens to be
%smaller than this minimum, it will be increased to the minimum.
if nargin<4 || isempty(minRadius)
    error('--plusTipTracker: min search radius input missing')
else
    parameters.minSearchRadius = minRadius;
end
%maximum allowed search radius. if a feature's calculated search radius is
%larger than this maximum, it will be reduced to this maximum.
if nargin<5 || isempty(maxRadius)
    error('--plusTipTracker: max search radius input missing')
else
    parameters.maxSearchRadius = maxRadius;
end


parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.
parameters.linearMotion = 1; %use linear motion Kalman filter.
parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
parameters.nnWindow = gapCloseParam.timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).

parameters.kalmanInitParam.searchRadiusFirstIteration = 20; %Kalman filter initialization parameters.
%parameters.kalmanInitParam = [];

costMatrices(1).parameters = parameters;
clear parameters

%% cost matrix for gap closing

%function name
costMatrices(2).funcName = 'costMatLinearMotionCloseGaps_EB3';

if nargin<6 || isempty(maxFAngle)
    parameters.maxFAngle = 30;
else
    parameters.maxFAngle = maxFAngle;
end

if nargin<7 || isempty(maxBDist)
    parameters.backVelMultFactor = 1.5;
else
    parameters.backVelMultFactor = maxBDist;
end

costMatrices(2).parameters = parameters;
clear parameters

%% Kalman filter function names

kalmanFunctions.reserveMem = 'kalmanResMemLM';
kalmanFunctions.initialize = 'kalmanInitLinearMotion_EB3';
kalmanFunctions.calcGain = 'kalmanGainLinearMotion_EB3';
kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';

%% additional input

%saveResults
trackDir = [runInfo.anDir filesep 'track'];
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

