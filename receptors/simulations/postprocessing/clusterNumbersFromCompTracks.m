function clusterStats = clusterNumbersFromCompTracks(compTracks,observeSideLen,probDim)
%clusterNumbersFromCompTracks calculates number, fraction and density of each cluster size from compTracks.
%
%   Input:
%       1) compTracks: The tracks structure array in default format.
%       2) observeSideLen:  Simulation/image side length values,
%                           which can be a single value or a value per side.
%       3) probDim: Problem dimension (default: 2).
%
%   Output:
%       clusterStats: a struct with the following fields
%           a) clusterCount (max cluster x num iters x num sims): number of 
%              clusters per size per iteration per simulation
%           b) clusterFrac (max cluster x num iters x num sims): fraction
%              of clusters per size per iteration per simulation
%           c) clusterDensity (max cluster x num iters x num sims): density
%              of clusters per size per iteration per simulation
%           d) receptorDensity (num sims x num iters): density of receptors
%              per iteration per simulation
%           e) largestClusterSize: the largest cluster size from all sims
%           f) observeSideLen: the (input) value used for the calculations
%           g) probDim: the (input) value used for the calculations
%           h) simArea: calculated area of the simulation space used when
%              determining densities.
%
%   Khuloud Jaqaman, May 2015. Based on calcClusterStatsFromCompTracks.
%

clusterStats = [];

if (nargin == 1)
    observeSideLen = [];
    probDim = [];
    simArea = NaN;
elseif (nargin == 2)
    %Assign default value when not provided
    probDim = 2;
end

numSideLenVals = length(observeSideLen);

%Determine area if possible
if (numSideLenVals == 1)
    simArea = observeSideLen ^ probDim;
elseif (numSideLenVals > 1)
    simArea = prod(observeSideLen);
end

%Get number of tracks and number of frames/iterations
numTracks = length(compTracks);
seqOfEvents = vertcat(compTracks.seqOfEvents);
numIters = max(seqOfEvents(:,1));

%reserve memory using a very large cluster size to be safe
maxClustSize = 1000;
clusterCount = zeros(maxClustSize,numIters);

%go over each track and count clusters at each time point
for iTrack = 1 : numTracks
    
    %get track's aggregation state
    aggregStateTrack = compTracks(iTrack).aggregState;
    
    %count clusters
    clustCountTmp = hist(aggregStateTrack,0:maxClustSize);
    
    %get the time span that this track covers
    colStart = compTracks(iTrack).seqOfEvents(1,1);
    colEnd = compTracks(iTrack).seqOfEvents(end,1);
    
    %add track statistics in overall matrix
    clusterCount(:,colStart:colEnd) = clusterCount(:,colStart:colEnd) + ...
        clustCountTmp(2:end,:);
    
end

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
clusterStats.observeSideLen = observeSideLen;
clusterStats.probDim = probDim;
clusterStats.simArea = simArea;
