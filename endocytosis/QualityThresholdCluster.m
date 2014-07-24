function [clusteredParticles,clusterCentroids] = QualityThresholdCluster(particlePositions,minClusterSize,clusterRadius);

% QualityThresholdCluster clusters the particles defined by their positions
% in input particlePositions into clusters of minimum size minClusterSize
% and maximum inter-particle distance clusterRadius
%
% INPUT
%           particlePositions = position of particles in two columns; first
%               column represents x-coordinate and second y-coordinate
%           minClusterSize = minimum number of particles that can be in a
%               cluster; defines minimum size of cluster
%           clusterRadius = maximum distance in between all particles
%               in cluster
% OUTPUT
%           clusterParticles = same matrix as particlePositions input but
%               with an added column that contains an id number for each
%               particle, so that a particle at (x1,y1) belongs to cluster
%               that can be identified by the number in the third column of
%               row 1.
%           clusterCentroids = position of cluster centroids in two
%               columns; first column is x position and second column is y
%               position; 
%
% REMARKS   The clustering is an adaptation of a gene clustering algorithm
% termed Quality Threshold Clustering. It assigns to all particles nearest
% neighbors until the maximum inter-particle distance is surpased. It then
% chooses the most unique. If there is no single most unique cluster it
% chooses amongst the most unique the one with the minimum average particle
% to centroid distance. If a single cluster can still not be identified it
% chooses amongst the pits with minimum particle to centroid distance the
% cluster with the minimum particle to cluster distance standard deviation.
%
% Uses: 
%
% Daniel Nunez, August 6, 2008

%add a column to particle position matrix which will hold an id number for
%the cluster to which particle is assgined
clusteredParticles = [particlePositions zeros(size(particlePositions,1),1)];
%calculate pdist of particles to be clustered
distMat = squareform(pdist(particlePositions));

%used to keep track of cluster ID
clusterNumber = 0;
%used to recod cluster centroids
clusterCentroids = [];

%repeat clustering until no cluster of min cluster size or larger can be
%found
maxClusterSize = minClusterSize;
tic
while maxClusterSize > minClusterSize - 1

    %initiate cluster descriptors that will be used later to pick the
    %best cluster (since particles belonging to a cluster are erased after
    %cluster is assigned and clusters are reassigned, these
    %decriptors need to be erased and recalculated everytime a cluster
    %is chosen)
    clusterSTD = [];
    clusterCentroid = [];
    clusterSize = [];
    clusterAvg = [];
    clusterList = [];

    %make cluster for each particle
    for iparticle = 1:length(distMat)
        if ~all(isnan(distMat(iparticle,:)))
            %reset cluster specific values
            clusterDiameter = 0;
            clusterParticles = [];
            distances = distMat(iparticle,:);
            %while cluster diamter is less than twice the max cluster radius specified
            while isempty(clusterDiameter) || clusterDiameter < 2*clusterRadius && ~all(isnan(distances))
                %find nearest neighbor
                findNearestNeigh = find(distances == min(distances));
                %add nearest neighbor to cluster
                clusterParticles = [clusterParticles findNearestNeigh];
                %calculate cluster diameter
                clusterDiameter = unique(max(pdist(particlePositions(clusterParticles,:))));
                
                if clusterDiameter > 2*clusterRadius
                %take out last value since this one puts cluster over limit radius
                clusterParticles(end) = [];
                end
                
                %erase neighbor from distance vector
                distances(findNearestNeigh) = nan;
            end %while clustersize is smaller than twice cluster radius

            %calculate size of cluster in terms of number of pits and add
            %to cluster list
            clusterSize(iparticle) = length(clusterParticles);
            %save clustered particle index numbers for future reference
            clusterList(iparticle).particles = clusterParticles;
        end %of if particle has not been assigned
    end %for each particle

    %find clusters with the greatest number of pits
    maxClusterSize = unique(max(clusterSize));
    findCluster = find(clusterSize == maxClusterSize);

    %sort more complicated clusters
    if maxClusterSize > minClusterSize - 1
        %if multiple clusters find overlap of particels in cluster, std from
        %centroid of cluster, and particle/centroid ditances for each cluster
        for ifind = 1:length(findCluster)
            %get index for pits in cluster
            particleID = clusterList(findCluster(ifind)).particles;
            %calculate centroid
            clusterCentroid(ifind,:) = sum(particlePositions(particleID,1:2),1)/length(particleID);
            %calculate std of distance from points to centroid
            clusterSTD(ifind) = std(sqrt((particlePositions(particleID,1)-clusterCentroid(ifind,1)).^2 + (particlePositions(particleID,2) - clusterCentroid(ifind,2)).^2));
            %find avg distance to centroid
            clusterAvg(ifind) = nanmean(sqrt((particlePositions(particleID,1)-clusterCentroid(ifind,1)).^2 + (particlePositions(particleID,2) - clusterCentroid(ifind,2)).^2));
        end %for each found cluster
        %out of these find the ones with the closest pits
        findMinAvg = find(round(10^4*clusterAvg) == min(round(10^4*clusterAvg)));
        %out of these find the ones with the min std of cluster/pit
        %distances
        findMinSTD = find(round(10^4*clusterSTD(findMinAvg)) == min(round(10^4*clusterSTD(findMinAvg))));
        %pick out final chosen cluster
        findCluster = findCluster(findMinAvg((findMinSTD)));
        %just in case there are still two clusters that cannot be
        %differentiated
        %if all remaining clusters are made up of same particles then
        %cluster is valid
        if length(findCluster)>1 & isempty(setdiff(clusterList(findCluster(1)).particles,[clusterList(findCluster(2:end)).particles]));
            findCluster = findCluster(1);
            findMinSTD = findMinSTD(1);
            %if one cluster has a particle not present in the others then this issue needs to be revisited
        elseif length(findCluster)>1 & ~isempty(setdiff(alteredPos(findNearNeigh(findCluster(1)).pits),alteredPos(pitSetList)))
            error('can not differetiate clusters which are different')
        end

        %assign cluster ID value to particles in cluster
        clusterNumber = clusterNumber + 1;
        clusteredParticles(clusterList(findCluster).particles,3) = clusterNumber;
        %calculate and record cluster centroids
        clusterCentroids = [clusterCentroids; clusterCentroid(findMinAvg(findMinSTD),:)];
        %erase particles in chosen cluster from further consideration in
        %clustering
        distMat(:,[clusterList(findCluster).particles]) = nan(size(distMat,1),length(clusterList(findCluster).particles));
        distMat([clusterList(findCluster).particles],:) = nan(length(clusterList(findCluster).particles),size(distMat,2));
    end %of if maxClusterSize is large enough

end %of while particles left to assign
toc
end %of function
