function plusTipCometTracker(projData,timeWindow,minTrackLen,minRadius,maxRadius,maxFAngle,maxBAngle,maxShrinkFactor,fluctRad,timeRange,diagnostics)
% plusTipCometTracker is the main tracking function


% get projData in correct format
if nargin<1 || isempty(projData)
    homeDir=pwd;
    projData.anDir=uigetdir(pwd,'Please select analysis directory');
    cd([projData.anDir filesep '..'])
    projData.imDir=uigetdir(pwd,'Please select image directory');
    cd(homeDir)
else
    % adjust for OS
    if ~isfield(projData,'imDir') || ~isfield(projData,'anDir')
        error('--plusTipCometTracker: first argument should be a structure with fields imDir and anDir');
    else
        [projData.anDir] = formatPath(projData.anDir);
        homeDir=pwd;
        cd(projData.anDir)
        [projData.imDir] = formatPath(projData.imDir);
        cd(homeDir)
    end
end

% load movieInfo (detection result)
featDir  = [projData.anDir filesep 'feat'];
if ~isdir(featDir)
    error('--plusTipCometTracker: feat directory missing')
else
    if exist([featDir filesep 'movieInfo.mat'],'file')
        load([featDir filesep 'movieInfo.mat'])
    else
        error('--plusTipCometTracker: movieInfo missing...')
    end
end


%% cost matrix for frame-to-frame linking

%function name
costMatrices(1).funcName = 'plusTipCostMatLinearMotionLink'; 

%minimum allowed search radius. the search radius is calculated on the spot
%in the code given a feature's motion parameters. if it happens to be
%smaller than this minimum, it will be increased to the minimum.
parameters.minSearchRadius = minRadius;

%maximum allowed search radius. if a feature's calculated search radius is
%larger than this maximum, it will be reduced to this maximum.
parameters.maxSearchRadius = maxRadius;

parameters.brownStdMult = 3; %multiplication factor to calculate search radius from standard deviation.
parameters.linearMotion = 1; %use linear motion Kalman filter.
parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
parameters.nnWindow = timeWindow; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).

parameters.kalmanInitParam.searchRadiusFirstIteration = 20; %Kalman filter initialization parameters.

% use movieInfo to check start/end frames for tracking
% then get rid of detected features in movieInfo outside this range
nFrames = length(movieInfo);
if nargin<10 || isempty(timeRange)
    startFrame=1;
    endFrame=nFrames;
elseif isequal(unique(size(timeRange)),[1 2])
    if timeRange(1)<=timeRange(2) && timeRange(2)<=nFrames
        startFrame = timeRange(1);
        endFrame = timeRange(2);
    else
        startFrame = 1;
        endFrame = nFrames;
    end
else
    error('--plusTipCometTracker: timeRange should be [startFrame endFrame] or [] for all frames')
end
temp=movieInfo;
clear movieInfo;
[movieInfo(1:nFrames,1).xCoord] = deal([]);
[movieInfo(1:nFrames,1).yCoord] = deal([]);
[movieInfo(1:nFrames,1).amp] = deal([]);
[movieInfo(1:nFrames,1).int] = deal([]);
[movieInfo(1:nFrames,1).ecc] = deal([]);
movieInfo(startFrame:endFrame,:)=temp(startFrame:endFrame,:);

% save the tracking frameRange
parameters.startFrame = startFrame;
parameters.endFrame = endFrame;

% run diagnostics on search radius range - get histogram of the linking distance 
% enter vector containing frames of interest for which you want to make the
% plot. keep in mind that 3 figures will be created for each frame, because
% the linking step is repeated forward, backward, and forward. used in
% plusTipCostMatLinearMotionLink.
if nargin<11 || isempty(diagnostics)
    parameters.diagnostics = [];
else
    parameters.diagnostics = round(mean([startFrame; endFrame]));
end

costMatrices(1).parameters = parameters;
clear parameters

%% cost matrix for gap closing

%function name
costMatrices(2).funcName = 'plusTipCostMatCloseGaps';

parameters.maxFAngle = maxFAngle;
parameters.maxBAngle = maxBAngle;
parameters.backVelMultFactor = maxShrinkFactor;
parameters.fluctRad = fluctRad;

costMatrices(2).parameters = parameters;
clear parameters

% maximum allowed time gap (in frames) between a track segment end and a
% track segment start that allows linking them. if timeWindow = n, then
% there are n-1 frames between the last detection event of the first track
% segment and the first detection event of the second track segment.
gapCloseParam.timeWindow = timeWindow;

% minimum length of track segments (in frames) from linking to be used in
% gap closing. 1 if everything, 2 if all tracks 2 or more frames long, etc.
gapCloseParam.minTrackLen = minTrackLen;

% 1 if merging and splitting are to be considered, 0 otherwise.
% even though EB comets might overlap in a given frame, they do not merge
% or split physically, so this should be zero.
gapCloseParam.mergeSplit = 0;

% run diagnostics on time window parameter - get histogram of all gap lifetimes
% use 1 to run and 0 to skip this step. used in trackCloseGapsKalman.
if nargin<11 || isempty(diagnostics)
    gapCloseParam.diagnostics = 0;
else
    gapCloseParam.diagnostics = 1;
end

%% Kalman filter function names

kalmanFunctions.reserveMem  = 'kalmanResMemLM';
kalmanFunctions.initialize  = 'plusTipKalmanInitLinearMotion';
kalmanFunctions.calcGain    = 'plusTipKalmanGainLinearMotion';
kalmanFunctions.timeReverse = 'kalmanReverseLinearMotion';

%% additional input

%saveResults
trackDir = [projData.anDir filesep 'track'];
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

[tracksFinal,kalmanInfoLink,errFlag,diagnosticTrackLinearity] = plusTipTrackCloseGapsKalman(movieInfo,...
    costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);

