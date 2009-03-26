function [clusterResults] = QTcluster(experiment,rest,force);

% QTcluster clusters initiations as defined by the rest vector in each
% movie in experiment and saves results under each movie directory.
%
% INPUT:   experiment=   data structure pointing to all the
%                       image/detection/tracking data for a given condition;
%                       for e.g. 10 cell movies, the structure should have 10
%                       entries; the data structure can be created with the
%                       function loadConditionData, and it needs to contain
%                       at least the field
%                       .source, which is the (path) location of the
%                       lifetime information folder
%                       .framerate, which is the movie framerate, which is
%                       necessary for the lifetime restriction
%           rest    =   restriction vector can have variable length;
%                       minimum length is five, where the entries are
%                       [stat da minfr minlft maxlft]
%                       optionally, the length can be extended to nine,
%                       where the additional entries are
%                       [... minint maxint minmot maxmot]
% OUTPUT
%           clusterResults.clusterResults = [xPos yPos clusterID lifetime startFrame]
%               the first column contains the x position, the second the y
%               position for all the initiations in a given movie
%               the third column contains an id number for each
%               particle, so that a particle at (x1,y1) belongs to cluster
%               that can be identified by the number in the third column of
%               row 1.
%           clusterResults.clusterCentroids = position of cluster centroids in two
%               columns; first column is x position and second column is y
%               position;
%           clusterResults.hotSpotRadius = radius used for clustering by
%               QualityThresholdCluster
%
% Uses:
%       determineHotSpotRadius
%       QualityThresholdCluster
%
% Daniel Nunez, updated March 11, 2009

%save old directory
oldDir = cd;
%set force variable
if exist('force','var') == 0 || isempty(force)
    force = 0;
end


%Fill in Missing Data
%movie length is required for
[experiment] = determineMovieLength(experiment);

%GET HOT SPOT RADIUS FROM DENSITY PLOTS
[hotSpotRadius] = determineHotSpotRadius(experiment,rest);


%%
%FOR EACH MOVIE
for iexp = 1:length(experiment)

    waitHandle = waitbar(iexp/length(experiment),['clustering in progress ' num2str(iexp) ' out of ' num2str(length(experiment))]);
    
    if exist([experiment(iexp).source filesep 'Cluster'],'dir') == 0 || force == 1 
    
    %Load Lifetime Information
    cd([experiment(iexp).source filesep 'LifetimeInfo'])
    load('lftInfo')
    % status matrix
    statMat = full(lftInfo.Mat_status);
    % lifetime matrix
    lftMat = full(lftInfo.Mat_lifetime);
    % x-coordinate matrix
    matX = full(lftInfo.Mat_xcoord);
    % y-coordinate matrix
    matY = full(lftInfo.Mat_ycoord);
    % disapp status matrix
    daMat = (lftInfo.Mat_disapp);
    % framerate
    framerate = experiment(iexp).framerate;

    %find all pits in movie that meet requirements specified by restriction
    %vector
    findPos = find((statMat==rest(1,1)) & (daMat==rest(1,2)) &...
        (lftMat>rest(1,3)) & (lftMat>round(rest(1,4)/framerate)) & (lftMat<round(rest(1,5)/framerate)));
    alteredPos = findPos;
    %find pits that are too far away from other pits to be clustered
    alteredDistMat = squareform(pdist([matX(alteredPos) matY(alteredPos)]));
    %make zeros into nans since min of this distance matrix will always be
    %zero otherwise, which is the distance from one pit to its own self
    alteredDistMat(alteredDistMat == 0) = nan;
    findLonelyPits = find(min(alteredDistMat,[],2) > 2*hotSpotRadius);
    %store pits that are too far from other pits as pits outside hotspots
    outsidePits = alteredPos(findLonelyPits);
    %erase these pits from alteredPos
    alteredPos(findLonelyPits) = [];

    particlePositions = [matX(alteredPos) matY(alteredPos)];
    
    %CLUSTER
    [clusteredParticles,clusterCentroids] = QualityThresholdCluster(particlePositions,2,hotSpotRadius);

    %add removed unclustered pits
    outsideParticles = [matX(outsidePits) matY(outsidePits) zeros(length(outsidePits),1)];
    clusteredParticles = [clusteredParticles; outsideParticles];

    %add lifetimes
    clusteredParticles(1:size(clusteredParticles,1),4) = lftMat([alteredPos;outsidePits])';
    %add start frame
    [dummy,startFrame] = ind2sub(size(lftMat),[alteredPos;outsidePits]);
    clusteredParticles(1:size(clusteredParticles,1),5) = startFrame';
    %make result structure for movie
    clusterResults.clusterResults = clusteredParticles;
    clusterResults.clusterCentroids = clusterCentroids;
    clusterResults.hotSpotRadius = hotSpotRadius;

    %save data onto folder under cell directory
    PATHNAME = experiment(iexp).source;
    cd(PATHNAME);
    %if ClusterData directory does not exist make it
    if ~(exist([PATHNAME filesep 'ClusterData'],'dir') == 7)
    mkdir(PATHNAME,'ClusterData')
    end
    %save cluster results under ClusterData folder; if cluster results
    %already exists, do not overwrite but rather add a number by which to
    %identiofy this particular clusterResults; functions that use
    %clusterResults will have to take the latest result, or have the user
    %specify a number to identy the result
    filePath = [PATHNAME filesep 'ClusterData' filesep 'clusterResults_rad' num2str(round(hotSpotRadius)) '_' num2str(round(rest(4)))...
        'to' num2str(round(rest(5))) 'lft'];
    securesave(filePath,'clusterResults');
    
    end %of if needs to cluster
    
    close(waitHandle)
    
end %of for each movie

%return to old directory
cd(oldDir)
end %of function
