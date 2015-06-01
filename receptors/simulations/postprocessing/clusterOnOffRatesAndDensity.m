function [rateOnPerClust,rateOffPerClust,densityPerClust,clustHistory,clustStats] = ...
    clusterOnOffRatesAndDensity(compTracksAggregState,infoSpaceTime)
%clusterOnOffRatesAndDensity calculates cluster on and off rates and densities
%
%   INPUT:   
%       compTracksAggregState: Compound tracks as output by
%                              aggregStateFromCompTracks_new. Contains
%                              tracks in both default and alternative
%                              formats.
%                            
%       infoSpaceTime: Structure with fields:
%           .probDim        : Problem dimensionality.
%           .areaSideLen    : Simulation/image side length values,
%                             which can be a single value or a value per
%                             side. In units of interest (e.g. um).
%           .timeStep       : Time between frames/time points. In units of
%                             interest (e.g. s).
%
%   OUTPUT:
%       rateOnPerClust    :  A 1D array of calculated on rates for clusters
%                            of size 1, 2, 3, etc. Cluster of size 1 gets
%                            NaN, but value kept for ease of reference to
%                            larger clusters. Units: per second per (#
%                            molecules/unit area). Area unit = square of
%                            areaSideLen unit.
%       rateOffPerClust   :  A 1D array of calculated off rates for clusters
%                            of size 1, 2, 3, etc. Cluster of size 1 gets
%                            NaN, but value kept for ease of reference to
%                            larger clusters. Units: per second.
%       densityPerClust   :  A 1D array of calculated density for clusters
%                            of size 1, 2, 3, etc. Units: # molecules/unit
%                            area. Area unit = square of areaSideLen unit.
%       clustHistory      : Output of
%                           clusterHistoryFromCompTracks_aggregState.
%       clustStats        : Output of clusterNumbersFromCompTracks.
%
%   Khuloud Jaqaman, May 2015

%% Input

%get space and time information
probDim = infoSpaceTime.probDim;
areaSideLen = infoSpaceTime.areaSideLen;
timeStep = infoSpaceTime.timeStep;

%% Calculation

%get cluster history
% *** 1st line (commented out) gets complete cluster history, including clusters
%that do not start or end in middle of time lapse; however it is very slow
% *** 2nd line (used for now) gets cluster history only for clusters that start
%AND end in midle of time lapse; this suffices for current analysis
% % [clustHistory,clustHistoryMerged] = clusterHistoryFromCompTracksComplete( ...
% %     compTracksAggregState.alternativeFormatTracks);
[clustHistory,clustHistoryMerged] = ...
    clusterHistoryFromCompTracks_aggregState(compTracksAggregState.defaultFormatTracks);

%get cluster densities
clustStats = clusterNumbersFromCompTracks(compTracksAggregState.defaultFormatTracks,areaSideLen,probDim);
densityPerClust = mean(clustStats.densityPerClust,2);

%get maximum cluster size
maxClusterSize = max(clustHistoryMerged(:,2));
rateOffPerClust = NaN(maxClusterSize,1);
rateOnPerClust = NaN(maxClusterSize,1);

%go over each cluster size > 1 and calculate off rate
for iSize = 2 : maxClusterSize
    
    %get lifetimes of clusters of current size, and whether they ended by
    %association or dissociation
    %only look at clusters with known start and end time
    indxClust = find(clustHistoryMerged(:,2)==iSize&~isnan(clustHistoryMerged(:,5)));
    clustLft = clustHistoryMerged(indxClust,5);
    clustEndType = clustHistoryMerged(indxClust,6);
    
    %calculate dissociation rate
    rateOffPerClust(iSize) = ...
        (length(find(clustEndType==1))/length(clustEndType)) / ...
        (mean(clustLft)*timeStep);
    
end

%calculate on rates from off rates and densities (assumes steady state)
rateOnPerClust(2:end) = rateOffPerClust(2:end) .* densityPerClust(2:maxClusterSize) ./ ...
    ( densityPerClust(1:maxClusterSize-1) * densityPerClust(1) );

%% ~~~ the end ~~~