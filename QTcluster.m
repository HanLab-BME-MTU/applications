function [experiment] = QTcluster(experiment,rest,force);

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

if nargin < 2 || isempty(rest)
    rest = [1 1 4 1 300];
end
%set force variable
if nargin < 3 || isempty(force)
    force = 0;
end


%GET HOT SPOT RADIUS FROM DENSITY PLOTS
[hotSpotRadius] = determineHotSpotRadius(experiment,rest);
%hotSpotRadius = 2.52;
%warning('hotSpotRadius preSet')
%%
%FOR EACH MOVIE
for iexp = 1:length(experiment)
        
    if exist([experiment(iexp).source filesep 'ClusterData'],'dir') == 0 || force == 1
        
        % framerate
        framerate = experiment(iexp).framerate;
        
        load([experiment(iexp).source filesep 'Tracking' filesep 'trackAnalysis.mat'])
        
        %track.status == 1 means track is complete
        tracks = tracks([tracks.status] == 1 & arrayfun(@(x)all(x.gapStatus~=5),tracks) == 1 &...
            [tracks.lifetime_s] > rest(1,3)*framerate & ...
            [tracks.lifetime_s] > rest(1,4) & [tracks.lifetime_s] < rest(1,5));
        
        particlePositionsX = arrayfun(@(t) t.x(1),tracks)';
        particlePositionsY = arrayfun(@(t) t.y(1),tracks)';
        
        %find all pits in movie that meet requirements specified by restriction
        %vector
        %findPos = find((statMat==rest(1,1)) & (daMat==rest(1,2)) &...
        %(lftMat>rest(1,3)) & (lftMat>round(rest(1,4)/framerate)) & (lftMat<round(rest(1,5)/framerate)));
        
        %find pits that are too far away from other pits to be clustered
        alteredDistMat = squareform(pdist([particlePositionsX particlePositionsY]));
        %make zeros into nans since min of this distance matrix will always be
        %zero otherwise, which is the distance from one pit to its own self
        alteredDistMat(alteredDistMat == 0) = nan;
        findLonelyPits = find(min(alteredDistMat,[],2) > 2*hotSpotRadius);
        
        lonelyParticles = arrayfun(@(x) [particlePositionsX(x) particlePositionsY(x)...
            0 tracks(x).lifetime_s/framerate tracks(x).start],findLonelyPits,'UniformOutput',false);
        lonelyParticles = cell2mat(lonelyParticles);
        
        tracks(findLonelyPits) = [];
        particlePositionsX(findLonelyPits) = [];
        particlePositionsY(findLonelyPits) = [];
        
        
        %CLUSTER
        [clusteredParticles,clusterCentroids] = QualityThresholdCluster([particlePositionsX particlePositionsY],2,hotSpotRadius);
        
        %add lifetimes and startframe
        lftAndFrames = arrayfun(@(x) [x.lifetime_s/framerate x.start],tracks','UniformOutput',false);
        clusteredParticles(1:size(clusteredParticles,1),4:5) =  cell2mat(lftAndFrames);
        
        %add lonelyPits
        clusteredParticles = [clusteredParticles; lonelyParticles];
        
        %make result structure for movie
        clusterResults.clusterResults = clusteredParticles;
        clusterResults.clusterCentroids = clusterCentroids;
        clusterResults.hotSpotRadius = hotSpotRadius;
        
        %save data onto folder under cell directory
        PATHNAME = experiment(iexp).source;

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
            'to' num2str(round(rest(5))) 'lft_' datestr(now, 'yyyymmdd')];
        secureSave(filePath,'clusterResults');
        
        experiment(iexp).clusterResults = clusterResults;
        
    end %of if needs to cluster
    
end %of for each movie

end %of function
