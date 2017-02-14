function clusterStatsLastFrame = clusterDensityLastFrame(compTracksAggregState,infoSpaceTime)
%CLUSTERDENSITYLASTFRAME calculates cluster densities in the last frame
%
%   SYNPOSIS: clusterStatsLastFrame = clusterDensityLastFrame(movieInfoLastFrame,infoIntensitySpace)
%
%   INPUT:   
%
%       compTracksAggregState: in the static case is a list of detected 
%                              features in the last frame, as output by 
%                              genMovieInfoFromTracksTracksSparse
%
%
%       infoSpaceTime: Structure with fields:
%                   
%                   .intensityInfo: Row vector with unit intensity mean 
%                    and standard deviation (e.g. the intensity of a single
%                    fluorophore labeling a single receptor).
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
intensityInfo = infoSpaceTime.intensityInfo(1);
probDim = infoSpaceTime.probDim;
areaSideLen = infoSpaceTime.areaSideLen;


%% Calculate the density

% load the intensity information for all segments present in the last frame

intensityVector=compTracksAggregState.amp(:,1);


%divide intensity by unit intensity to get cluster size
clustSizeVec = round(intensityVector/intensityInfo);
clustSizeVec(clustSizeVec==0) = 1;


% calculate the maximum cluster size

maxClustSize=max(clustSizeVec);

%use hist function to count number of clusters of different sizes
clusterCount = hist(clustSizeVec,1:maxClustSize);

% %reserve memory, the number of row is equivalent to cluster size and for
% %each cluster size the number of elements in that
% 
% clusterCount = zeros(maxClustSize,1);
% 
% %For each segment in movieInfoLastFrame calculate the cluster size
% 
% for segIndex=1:length(intensityVector)
% % calculate the round value of intensity, 1 cluster size 1, 2 cluster size
% % 2 and so on
% 
% intensityValue=round(intensityVector(segIndex)/intensityInfo);
% if intensityValue~=0
% clusterCount(intensityValue) = ...
% clusterCount(intensityValue) + 1;
% end
% end


% to have the same outputs as in the function clusterNumbersFromCompTracks

%Determine area to be used in the calculation of density per cluster
numSideLenVals = length(areaSideLen);
if (numSideLenVals == 1)
    simArea = areaSideLen ^ probDim;
else
    simArea = prod(areaSideLen);
end

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


