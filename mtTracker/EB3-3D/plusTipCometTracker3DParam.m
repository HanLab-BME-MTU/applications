function [costMatrices,gapCloseParam,kalmanFunctions,probDim]=plusTipCometTracker3DParam(MD)

%% Original tunable 
%% Tracking parameters
% Segments presenting a lifetime below <minTrackLength> are discarded. 
minTrackLen=3; % In frames. 3
 
% Lower and upper bound for the radius of the sphere considered to link 
% a detection one frame to the next
minRadius=2; % In pixels.
maxRadius=8; % In pixels. % 8
searchRadiusMult=3;
searchRadiusFirstIteration=10; %10

% In the context of pause or shrinkage detection, <maxFAngle> is the
% maximum angle between the estimated speed vectors at the end
% of a segment and a possible link between two segments. Also the maximum angle
% between the speed estimated the beginning and the end of a segment.
maxFAngle=10; % In degrees.

% In the context of shrinkage detection: <maxBAngle>  defines the area
% considered  around the segment to detect shrinkage event. The angle only
% defines the distance orthoganal distance between a point and track (as if
% the track was a straight line). See Figure 5B in Applegate and al 2012
% for an illustration.
maxBAngle=10; % In degrees.

% Maximum skrinkage factor (multiply with the maximum growth rate) 
maxShrinkFactor=1.5;

timeWindow=2;

% The fluctuation radius <fluctRad> models unexpected fluctuations during
% shrinkage or pauses. The search volume is defined bye area described by
% the parameter above and dilatted by the fluctuation radius.
fluctRad=1.; % In pixels.

% Break non linear tracks into multiple tracks.
breakNonLinearTracks=false;
diagnostics=[];

timeRange=[1 MD.nFrames_];

timeWindowNN=timeWindow;


%% cost matrix for frame-to-frame linking

%function name
costMatrices(1).funcName = func2str(@plusTipCostMatLinearMotionLink); 

%minimum allowed search radius. the search radius is calculated on the spot
%in the code given a feature's motion parameters. if it happens to be
%smaller than this minimum, it will be increased to the minimum.
parameters.minSearchRadius = minRadius;

%maximum allowed search radius. if a feature's calculated search radius is
%larger than this maximum, it will be reduced to this maximum.
parameters.maxSearchRadius = maxRadius;

parameters.brownStdMult = searchRadiusMult; %multiplication factor to calculate search radius from standard deviation.
parameters.linearMotion = 1; %use linear motion Kalman filter.
parameters.useLocalDensity = 1; %1 if you want to expand the search radius of isolated features in the linking (initial tracking) step.
parameters.nnWindow = timeWindowNN; %number of frames before the current one where you want to look to see a feature's nearest neighbor in order to decide how isolated it is (in the initial linking step).

parameters.kalmanInitParam.searchRadiusFirstIteration = searchRadiusFirstIteration; %Kalman filter initialization parameters.

% use movieInfo to check start/end frames for tracking
% then get rid of detected features in movieInfo outside this range
nFrames = MD.nFrames_;
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

% Initialize a new movieInfo structure array with the old movieInfo fields 
% movieInfo fields may vary depending of the detection method 
% oldMovieInfo=movieInfo;
% clear movieInfo;
% movieFields = fieldnames(oldMovieInfo)';
% emptyFields=[movieFields; cell(size(movieFields))];
% movieInfo(nFrames,1) = struct(emptyFields{:});
% movieInfo(startFrame:endFrame,:)=oldMovieInfo(startFrame:endFrame,:);
% clear oldMovieInfo

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
costMatrices(2).funcName = func2str(@plusTipCostMatCloseGaps);

parameters.maxFAngle = maxFAngle;
parameters.maxBAngle = maxBAngle;
parameters.backVelMultFactor = maxShrinkFactor;
parameters.fluctRad = fluctRad;
parameters.breakNonLinearTracks = breakNonLinearTracks;

    
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
    gapCloseParam.diagnostics = diagnostics;
end

%% Kalman filter function names

kalmanFunctions.reserveMem  = func2str(@kalmanResMemLM);
kalmanFunctions.initialize  = func2str(@plusTipKalmanInitLinearMotion);
kalmanFunctions.calcGain    = func2str(@plusTipKalmanGainLinearMotion);
kalmanFunctions.timeReverse = func2str(@kalmanReverseLinearMotion);

%% additional input

%verbose
verbose = 1;

%problem dimension
probDim = 3;






