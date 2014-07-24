function assignPointsToModels2(obj)
% function assignPointsToModels2(obj)
% SYNOPSIS:
% Assigns every point in the data set to the likeliest model.
%
% REQUIRED INPUTS:         
% 
% OPTIONAL INPUTS:
%
% NEEDED PROPERTIES: 
% - obj.data.points
% - obj.data.nPoints
% - obj.data.clusters
% - obj.data.nullCluster
% - obj.data.neighbors
% - obj.data.parents
% - obj.data.error
% - obj.data.modelBezCP
% - obj.data.modelType
% - obj.data.modelIsOutOfDate
% - obj.data.modelVar
%
% MODIFIED PROPERTIES: 
% - obj.data.clusters
% - obj.data.nullCluster
% - obj.data.modelIsOutOfDate
%
% OUTPUTS:
%
% Pascal Bérard, February 2012

% Compute the classification mixture component weigths
nPointsInCluster = cellfun(@numel,obj.data.clusters);
compWeights = nPointsInCluster/obj.data.nPoints;
nullCompWeight = numel(obj.data.nullCluster)/obj.data.nPoints;

% Create cluster IDs
clusterIDs = cellfun(@sort,obj.data.clusters,'UniformOutput',false);

% Backup the current clusters
clustersOld = obj.data.clusters;

% Remove the neighbors which are unclustered points
neighbors = cellfun(@(a) setdiff(a,obj.data.nullCluster),obj.data.neighbors,'UniformOutput',false);

% Get the parent cluster of the neighbor points
p = obj.data.parents;
parents = cellfun(@(a) p(a(:)),neighbors,'UniformOutput',false);

% Remove double entries
parents = cellfun(@unique,parents,'UniformOutput',false);

% Compute an estimate of the data volume
dim = any(obj.data.points ~= 0,1); % Dimensionality of the points
vol = prod(diff(quantile(obj.data.points(:,dim),[0.05,0.95],1))./0.9); %*400/500;
% vol = prod(max(obj.data.points(:,dim),[],1)-min(obj.data.points(:,dim),[],1));

% Determine the likelihood threshold
logLThreshold = log(1/vol*nullCompWeight);

likeliestParent = cell(obj.data.nPoints,1);

obj_data_modelBezCP = obj.data.modelBezCP;
obj_data_error = obj.data.error;
obj_data_modelType = obj.data.modelType;
obj_data_points = obj.data.points;
obj_data_modelVar = obj.data.modelVar;
obj_data_modelLength = obj.data.modelLength;

parfor i=1:obj.data.nPoints
    % Compute the likelihood
    dist = arrayfun(@(b) distancePointBezier(obj_data_modelBezCP{b}./repmat(obj_data_error(i,:),obj_data_modelType(b)+1,1), ...
        obj_data_points(i,:)./obj_data_error(i,:)), ...
        parents{i});
    modelLength = obj_data_modelLength(parents{i});
    sigmaComponent = sqrt(obj_data_modelVar(parents{i})/2);
    compWeightsNeighborModels = compWeights(parents{i});
    logL = arrayfun(@(a,c,d) Processor.logLikelihoodModelAlongAway2D([a 0 0],1./obj_data_error(i,:),c,d),dist.*obj_data_error(i,1),sigmaComponent,modelLength);
    weightedLogL = log(compWeightsNeighborModels)+logL;
    
    % Find the indices of the models with a bigger likelihood than the threshold likelihood
    idx = find(weightedLogL>logLThreshold);
    
    % Remove these models from the likelihood and parents arrays
    weightedLogL = weightedLogL(idx);
    parents{i} = parents{i}(idx);
    
    % Find the most likeliest model
    [~,idx] = max(weightedLogL);
    
    % Get the cluster index for the likeliest model
    likeliestParent{i} = parents{i}(idx);
end

% Add the point to the unclustered points if no closest model exists
emptyLikeliestParent = cellfun(@isempty,likeliestParent);
pointsIdx = 1:obj.data.nPoints;

% Rebuild cluster list
obj.data.clusters = cell(obj.data.nClusters,1);

for p=pointsIdx(~emptyLikeliestParent)
    obj.data.clusters(likeliestParent{p}) = {[obj.data.clusters{likeliestParent{p}} pointsIdx(p)]};
end

% Create new cluster IDs
newClusterIDs = cellfun(@sort,obj.data.clusters,'UniformOutput',false);

% Find clusters whose members changed
obj.data.modelIsOutOfDate = ~cellfun(@isequal,clusterIDs,newClusterIDs);

% Restore point order for the clusters that didn't change
obj.data.clusters(~obj.data.modelIsOutOfDate) = clustersOld(~obj.data.modelIsOutOfDate);

% Build null cluster
obj.data.nullCluster = permute(pointsIdx(emptyLikeliestParent),[2 1]);

disp('Process: All points have been assigned!');

end

