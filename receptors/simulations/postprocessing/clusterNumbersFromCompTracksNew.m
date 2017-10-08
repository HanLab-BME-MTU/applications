function clusterStats = clusterNumbersFromCompTracksNew(clustHistoryMerged,infoSpaceTime)
%CLUSTERNUMBERSFROMCOMPTRACKS calculates number, fraction and density of each cluster size from compTracks.
%
%   SYNPOSIS: clusterStats = clusterNumbersFromCompTracks(compTracks,infoSpaceTime)
%
%   Input:
%       clustHistoryMerged : Same information as in clustHistoryAll but with
%                         all cells merged into one 2D array, i.e.
%                         individual track information is lost.
%
%       infoSpaceTime: Structure with fields:
%           .probDim        : Problem dimensionality.
%           .areaSideLen    : Simulation/image side length values,
%                             which can be a single value or a value per
%                             side. In units of interest (e.g. um).
%           .timeStep       : Time between frames/time points. In units of
%                             interest (e.g. s).
%           .sampleStep     : Sampling time step, in same units as
%                             timeStep. Mostly relevant for
%                             simulated data where simulation time step
%                             might be 0.01 s but sampling time step of
%                             interest is e.g. 0.1 s.
%                             Optional. If not input, then sampleStep =
%                             timeStep.
%           .firstLastTP    : Row vector of first and last time points to
%                             use for calculating rates and densities. In
%                             same units as timeStep. 
%                             If only one value is input, it is taken as
%                             the last time point.
%                             If no value is input, then all time points
%                             are used.
%
%   Output:
%       clusterStats: a struct with the following fields
%           a) clusterCount (max cluster x num time points): number of 
%              clusters per size per time point.
%           b) clusterFrac (max cluster x num time points): fraction
%              of clusters per size per time point.
%           c) clusterDensity (max cluster x num time points): density
%              of clusters per size per time point.
%           d) receptorCount (1 x num time points): number of receptors
%              per time point.
%           e) receptorDensity (1 x num time points): density of receptors
%              per time point.
%           f) largestClusterSize: the largest cluster size.
%           g) infoSpaceTime: the (input) structure infoSpaceTime.
%
%   Khuloud Jaqaman, May 2015. Based on calcClusterStatsFromCompTracks.
%   Modified by Luciana de Oliveira, Oct 2016. To calculate the clustStats
%   directly from clustHistory.

%% Input

%get basic space and time information
probDim = infoSpaceTime.probDim;
areaSideLen = infoSpaceTime.areaSideLen;
timeStep = infoSpaceTime.timeStep;

%get sampling and cropping information
if isfield(infoSpaceTime,'sampleStep')
    sampleStep = infoSpaceTime.sampleStep;
else
    sampleStep = timeStep;
end
convStep = round(sampleStep/timeStep);
if isfield(infoSpaceTime,'firstLastTP')
    firstLastTP = infoSpaceTime.firstLastTP;
    if length(firstLastTP)==1
        firstLastTP = [0 firstLastTP];
    end
else
    firstLastTP = [];
end
firstLastFrame =1+  round(firstLastTP/sampleStep);


%Determine area
numSideLenVals = length(areaSideLen);
if (numSideLenVals == 1)
    simArea = areaSideLen ^ probDim;
else
    simArea = prod(areaSideLen);
end

%% Calculations

%Get number of tracks and number of frames/iterations
numOfEvents = length(clustHistoryMerged(:,1));
numIters = max(clustHistoryMerged(:,4));

%define default first and last frame if not input
if isempty(firstLastTP)
    firstLastFrame = [1 round(numIters/convStep)];
end

%reserve memory using a very large cluster size to be safe
maxClustSize = 1000;
clusterCount = zeros(maxClustSize,numIters);

%go over each track and count clusters at each time point
for iTrack =  1 : numOfEvents
    
    %get events from clustHistory
    eventTrack = clustHistoryMerged(iTrack,:);
    
    %get the time span that this track covers

    colStart = eventTrack(3);
    colEnd = eventTrack(4);
    
    % LRO, 2016/11/02: this is necessary to avoid the problem that in clust
    % history the same frame that end one segment can be the start of
    % another one
    
    if colEnd~=numIters
     colEnd = eventTrack(4)-1;
    end
    
    %count clusters and add to overall matrix
        
        clusterCount(eventTrack(2),colStart:colEnd) = ...
        clusterCount(eventTrack(2),colStart:colEnd) + 1;
           
end

%sub-sample time step and crop time as requested
clusterCount = clusterCount(:,1:convStep:end);
clusterCount = clusterCount(:,firstLastFrame(1):firstLastFrame(2));
numIters = size(clusterCount,2);

%remove rows after largest cluster size
numClustTot = sum(clusterCount,2);
largestClustSize = find(numClustTot~=0,1,'last');
clusterCount = clusterCount(1:largestClustSize,:);

%calculate cluster fractions at each iteration
clusterFrac = clusterCount ./ repmat(sum(clusterCount),largestClustSize,1);

%calculate cluster density at each iteration
clusterDensity = clusterCount / simArea;

%calculate receptor density at each iteraction
receptorCount = sum( clusterCount .* repmat((1:largestClustSize)',1,numIters) );
receptorDensity = receptorCount / simArea;

%Save values for return
clusterStats.clusterCount = clusterCount;
clusterStats.clusterFrac = clusterFrac;
clusterStats.clusterDensity = clusterDensity;
clusterStats.receptorCount = receptorCount;
clusterStats.receptorDensity = receptorDensity;

clusterStats.largestClustSize = largestClustSize;
clusterStats.infoSpaceTime = infoSpaceTime;

%% ~~~ the end ~~~

