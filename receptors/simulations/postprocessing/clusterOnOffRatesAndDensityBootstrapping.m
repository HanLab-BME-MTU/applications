function [rateOnPerClust,rateOffPerClust,densityPerClust,...
    numClustForRateCalc,clustStats] = ...
    clusterOnOffRatesAndDensityBootstrapping(clustHistoryMerged,infoSpaceTime)


%CLUSTERONOFFRATESANDDENSITY calculates cluster on and off rates and densities
% for each repetition of bootstrapping
%
%   SYNPOSIS: [rateOnPerClust,rateOffPerClust,densityPerClust,...
%     numClustForRateCalc,clustStats] = ...
% clusterOnOffRatesAndDensityBootstrapping(clustHistoryMerged,infoSpaceTime)
% 
% 
% INPUT:   
%      clustHistoryMerged: a 1D cell with rows = number of tracks in
%                         compTracks. 
%                         Each entry contains a clusterHistory table for a
%                         track in compTracks. In this version cluster 
%                         history is recorded for all clusters in the
%                         obervation time.
%                         clusterHistory is a 2D array with each row
%                         corresponding to an association or a dissociation
%                         event. The 9 colums give the following information:
%                         1) Track segment number.
%                         2) Cluster size.
%                         3) Start time (same units as input timeStep etc).
%                         4) End time (same units as input timeStep etc).
%                         5) Lifetime (same units as input timeStep etc).
%                         6) Event that started the cluster
%                           (1 = dissociation, 2 = association).
%                         7) Event that ended the cluster 
%                           (1 = dissociation, 2 = association).
%                         8) Resulting cluster size.
%                         9) Association flag - 1 indicates the segment
%                            and its partner are both listed and NaN
%                            indicates only the current segment is listed,
%                            i.e. the partner is not listed.
%                          All cells merged into one 2D array, i.e.
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
%
%   Modified by Luciana de Oliveira, Oct 2016. To calculate the clustStats
%   directly from clustHistory.

%% Input

%% Input

%get sampling information
timeStep = infoSpaceTime.timeStep;
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
firstLastFrame = 1 + firstLastTP/timeStep;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bootstrapping

%generates the random numbers:
numberOfRows=randi(size(clustHistoryMerged,1),size(clustHistoryMerged,1),1);

% bootstrap clustHystoryMerged
clustHistoryMerged=clustHistoryMerged(numberOfRows,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Densities
%get cluster densities
clustStats = clusterNumbersFromCompTracksNew(clustHistoryMerged,infoSpaceTime);

densityPerClust = mean(clustStats.clusterDensity,2);
clustCount = mean(clustStats.clusterCount,2);


%% RATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modification 2016/11/17(LRO): To calculate the rates it is needed to do the 
% following steps in clustHistoryMerged.

%Remove those events that are not changing

clustHistoryMergedChange=clustHistoryMerged;

 clustHistoryMergedChange(clustHistoryMergedChange(:,4)==max(clustHistoryMergedChange(:,4)),:) = [];
 clustHistoryMergedChange(clustHistoryMergedChange(:,3)==min(clustHistoryMergedChange(:,3)),:) = [];
    

 if ~isempty(firstLastFrame)
    clustHistoryMergedChange = clustHistoryMergedChange(clustHistoryMergedChange(:,3)>=firstLastFrame(1) ...
    & clustHistoryMergedChange(:,4)<=firstLastFrame(2),:);
 end
    
      
%remove clusters which start and end at exactly the same time point -
%these would be not detectable with the sub-sampling time step
    clustHistoryMergedChange = clustHistoryMergedChange(clustHistoryMergedChange(:,5)>0,:);
   
 %subsample time
    
   clustHistoryMergedChange(:,3:4) = ceil(clustHistoryMergedChange(:,3:4)/convStep);
    clustHistoryMergedChange(:,5) = clustHistoryMergedChange(:,4) - clustHistoryMergedChange(:,3);

%convert from iterations/frames to real time units
    clustHistoryMergedChange(:,3:5) = clustHistoryMergedChange(:,3:5) * sampleStep;    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%get maximum cluster size
maxClusterSize = min([max(clustHistoryMergedChange(:,2)) length(densityPerClust)]);
rateOffPerClust = NaN(maxClusterSize,1);
rateOnPerClust = NaN(maxClusterSize,1);
numClustForRateCalc = [NaN(maxClusterSize,1) clustCount(1:maxClusterSize)];

%go over each cluster size > 1 and calculate off rate
for iSize = 2 : maxClusterSize
    
    %get lifetimes of clusters of current size, and whether they ended by
    %association or dissociation
    %only look at clusters with known start and end time
    indxClust = find(clustHistoryMergedChange(:,2)==iSize&~isnan(clustHistoryMergedChange(:,5)));
    clustLft = clustHistoryMergedChange(indxClust,5);%cluster life time
    clustEndType = clustHistoryMergedChange(indxClust,7);
    
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


%% ~~~ the end ~~~