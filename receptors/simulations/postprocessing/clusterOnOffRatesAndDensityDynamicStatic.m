function [rateOnPerClust,rateOffPerClust,densityPerClust,...
    numClustForRateCalc,clustHistory,clustStats] = ...
    clusterOnOffRatesAndDensityDynamicStatic(compTracksAggregState,infoSpaceTime)
%CLUSTERONOFFRATESANDDENSITY calculates cluster on and off rates and densities
%
%   SYNPOSIS: [rateOnPerClust,rateOffPerClust,densityPerClust,...
%               numClustForRateCalc,clustHistory,clustStats] = ...
%               clusterOnOffRatesAndDensity(compTracksAggregState,infoSpaceTime)
%
%   INPUT:   
%       compTracksAggregState: Compound tracks as output by
%                              aggregStateFromCompTracks_new. Contains
%                              tracks in both default and alternative
%                              formats, including aggregation state.
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
% 
%            .systemState    : value saying if the analysis is for
%                              dynaminic or static data. Should be one of 
%                              the two values:
%                              1-dynamic;
%                              0-static;
%
%            .intensityInfo: Row vector with unit intensity mean 
%                    and standard deviation (e.g. the intensity of a single
%                    fluorophore labeling a single receptor).
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
%       numClustForRateCalc: First column indicates number of clusters of
%                            each size used to calculate off rate.
%                            Second column indicates mean number of
%                            clusters of each size per iteration,
%                            indirectly used to calculate on rate.
%       clustHistory      :  Output of clusterHistoryFromCompTracks_aggregState.
%       clustStats        :  Output of clusterNumbersFromCompTracks.
%
%   Khuloud Jaqaman, May 2015
%   Modified Luciana de Oliveira, February 2017. 
%% Input

%%%%% Modification LRO
%Test if it is dynamic or static data
systemState = infoSpaceTime.systemState;

if systemState==1
%% Dynamic Calculation

%get cluster history

% *** 1st line (commented out) gets complete cluster history, including clusters
%that do not start or end in middle of time lapse; however it is very slow
%ALSO, IT HAS A BUG FOR GETTING STARTING EVENT TYPES - FIX BEFORE USING
% % [clustHistory,clustHistoryMerged] = clusterHistoryFromCompTracksComplete( ...
% %     compTracksAggregState.alternativeFormatTracks);

% *** 2nd line (used for now) gets cluster history only for clusters that start
%AND end in middle of time lapse; this suffices for current analysis
[clustHistory,clustHistoryMerged] = ...
    clusterHistoryFromCompTracks_aggregState(compTracksAggregState.defaultFormatTracks,infoSpaceTime);

%get cluster densities
clustStats = clusterNumbersFromCompTracks(compTracksAggregState.defaultFormatTracks,infoSpaceTime);
densityPerClust = mean(clustStats.clusterDensity,2);
clustCount = mean(clustStats.clusterCount,2);

%get maximum cluster size
maxClusterSize = min([max(clustHistoryMerged(:,2)) length(densityPerClust)]);
rateOffPerClust = NaN(maxClusterSize,1);
rateOnPerClust = NaN(maxClusterSize,1);
numClustForRateCalc = [NaN(maxClusterSize,1) clustCount(1:maxClusterSize)];

%go over each cluster size > 1 and calculate off rate
for iSize = 2 : maxClusterSize
    
    %get lifetimes of clusters of current size, and whether they ended by
    %association or dissociation
    %only look at clusters with known start and end time
    indxClust = find(clustHistoryMerged(:,2)==iSize&~isnan(clustHistoryMerged(:,5)));
    clustLft = clustHistoryMerged(indxClust,5);%cluster life time
    clustEndType = clustHistoryMerged(indxClust,7);
    
    %calculate dissociation rate
    numClusters = length(clustEndType);
    rateOffPerClust(iSize) = ...
        (length(find(clustEndType==1))/numClusters) / mean(clustLft);
    
    %record number of clusters used for off rate calculation
    numClustForRateCalc(iSize,1) = numClusters;
    
end

%calculate on rates from off rates and densities (assumes steady state)
if maxClusterSize > 1
    rateOnPerClust(2:end) = rateOffPerClust(2:end) .* densityPerClust(2:maxClusterSize) ./ ...
        ( densityPerClust(1:maxClusterSize-1) * densityPerClust(1) );
end

%%%%%%%%%%%%%%% Modification LRO 2017/02/09

elseif systemState==0
%% static calculations
rateOnPerClust=[];
rateOffPerClust=[];
numClustForRateCalc=[];
clustHistory=[];
clustStats = clusterDensityLastFrame(compTracksAggregState,infoSpaceTime);
densityPerClust=clustStats.clusterDensity;
end

%% ~~~ the end ~~~