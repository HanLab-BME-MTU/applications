function filename = wyssFileTrackingSingle(data,name)
% takes a PointList from wyssFileImport and runs sparse kalman tracking on
% each point list.
%
%

%Set parameters for tracking

%% Creates parameter structures for tracking
%% general gap closing parameters
gapCloseParam.timeWindow = 3; %maximum allowed time gap (in frames) between a track segment end and a track segment start that allows linking them.
gapCloseParam.mergeSplit = 0; %1 if merging and splitting are to be considered, 2 if only merging is to be considered, 3 if only splitting is to be considered, 0 if no merging or splitting are to be considered.
gapCloseParam.minTrackLen = 1; %minimum length of track segments from linking to be used in gap closing.

%optional input:
gapCloseParam.diagnostics = 0; %1 to plot a histogram of gap lengths in the end; 0 or empty otherwise.

%% cost matrix for frame-to-frame linking

%function name
costMatrices(1).funcName = 'costMatStationaryLink';

%parameters
parameters.searchRadius = 2;

costMatrices(1).parameters = parameters;
clear parameters

%% cost matrix for gap closing

%function name
costMatrices(2).funcName = 'costMatStationaryCloseGaps';

%parameters
parameters.searchRadius = 2; 
parameters.gapPenalty = 2;

costMatrices(2).parameters = parameters;
clear parameters

%% Kalman filter function names

kalmanFunctions = [];

%% additional input
saveResults = 0;

%verbose
verbose = 1;

%problem dimension
probDim = 2;

%% Process anlyzed movies

    movieInfo = wyssFileConvertforTracking(data);
    
    %Apply tracking and Gap closing to localized data
    [tracksFinal,kalmanInfoLink,errFlag] = trackCloseGapsKalmanSparse(movieInfo,...
    costMatrices,gapCloseParam,kalmanFunctions,probDim,saveResults,verbose);

    features = zeros(size(movieInfo));

    clear movieInfo;
    
    %saves tracks into the current directory
    filename = [name(1:end-4),'_tracking.mat'];
    save(filename,'tracksFinal','features');    
    
    
end
