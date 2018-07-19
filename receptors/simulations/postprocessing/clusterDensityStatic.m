function clusterStatsLastFrame = clusterDensityStatic(detectionAggregState,infoSpaceTime)
%CLUSTERDENSITYLASTFRAME calculates cluster densities in the last frame
%
%   SYNPOSIS: clusterStatsLastFrame = clusterDensityLastFrame(movieInfoLastFrame,infoIntensitySpace)
%
%   INPUT:   
%
%       detectionAggregState: in the static case is a list of detected 
%                              features in the last frame, as output by 
%                              genMovieInfoFromTracksTracksSparse
%
%
%       infoSpaceTime: Structure with fields:                
%                  
%                    .probDim        : Problem dimensionality.
%
%                   .areaSideLen    : Simulation/image side length values,
%                             which can be a single value or a value per
%                             side. In units of interest (e.g. um).
%
%                   
%   OUTPUT:
% 
% 
%           clusterStats: a struct with the following fields (outputs as in
%           function clusterNumbersFromCompTracks)
% 
%           a) clusterCount (max clusterin in the last frame): number of 
%              clusters per size per time point.
%           b) clusterFrac (max cluster in the last frame): fraction
%              of clusters per size per time point.
%           c) clusterDensity (max cluster in the last frame): density
%              of clusters per size per time point.
%           d) receptorCount : number of receptors
%              in the last frame.
%           e) receptorDensity : density of receptors
%              in the last frame.

%% Input

%get intensity and space information
probDim = infoSpaceTime.probDim;
areaSideLen = infoSpaceTime.areaSideLen;

clusterCount=detectionAggregState;

% to have the same outputs as in the function clusterNumbersFromCompTracks

%Determine area to be used in the calculation of density per cluster
numSideLenVals = length(areaSideLen);
if (numSideLenVals == 1)
    simArea = areaSideLen ^ probDim;
else
    simArea = prod(areaSideLen);
end

maxClustSize=length(clusterCount);

%calculate cluster fractions 
clusterFrac = clusterCount ./ sum(clusterCount);

%calculate cluster density
clusterDensity = clusterCount / simArea;

%calculate receptor density per cluster size
receptorCountClust=zeros(length(clusterCount),1);
receptorDensityClust=zeros(length(clusterCount),1);

%number of receptors and clusters with size 1 is the same
receptorCountClust(1)=clusterCount(1);
receptorDensityClust(1)=clusterDensity(1);

for clusterIndex=2:maxClustSize
 receptorCountClust(clusterIndex) = clusterCount(clusterIndex).* clusterIndex;
 receptorDensityClust (clusterIndex)= receptorCountClust (clusterIndex)/ simArea;
end

% calculate the total number of receptors in the last frame

receptorCount=sum(receptorCountClust);
receptorDensity=sum(receptorDensityClust);


%Save values for return
clusterStatsLastFrame.clusterCount = clusterCount';
clusterStatsLastFrame.clusterFrac = clusterFrac';
clusterStatsLastFrame.clusterDensity = clusterDensity';
clusterStatsLastFrame.receptorCount = receptorCount;
clusterStatsLastFrame.receptorDensity = receptorDensity;

end


