function assignPointsToModels(obj,nSigmaThreshold)
% function assignPointsToModels(obj,nSigmaThreshold)
% SYNOPSIS:
% Assigns every point in the data set to the closest model.
%
% REQUIRED INPUTS:         
% - nSigmaThreshold
% If a point is further away than nSigmaThrehold times the error from its
% closest model, the point is classified as noise.
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
% Pascal Bérard, October 2011

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

% Determine the likelihood threshold for each point (The threshold is at nSigmaThreshold times the sqrt(component variance))
likelihoodThreshold = -0.5*nSigmaThreshold^2;

nearestParent = cell(obj.data.nPoints,1);

obj_data_modelBezCP = obj.data.modelBezCP;
obj_data_error = obj.data.error;
obj_data_modelType = obj.data.modelType;
obj_data_points = obj.data.points;
obj_data_modelVar = obj.data.modelVar;

parfor i=1:obj.data.nPoints
    likelihood = arrayfun(@(b)-(0.5*distancePointBezier(obj_data_modelBezCP{b}./repmat(obj_data_error(i,:),obj_data_modelType(b)+1,1), ...
        obj_data_points(i,:)./obj_data_error(i,:))^2/obj_data_modelVar(b)), ...
        parents{i});
    
    % Find the indices of the models with a bigger likelihood than the threshold likelihood
    idx = find(likelihood>likelihoodThreshold);
    
    % Remove these models from the likelihood and parents arrays
    likelihood = likelihood(idx);
    parents{i} = parents{i}(idx);
    
    % Find the most likeliest model
    [~,idx] = max(likelihood);
    
    % Get the cluster index for the likeliest model
    nearestParent{i} = parents{i}(idx);
end

% Add the point to the unclustered points if no closest model exists
emptyNearestParent = cellfun(@isempty,nearestParent);
pointsIdx = 1:obj.data.nPoints;

% Rebuild cluster list
obj.data.clusters = cell(obj.data.nClusters,1);

for p=pointsIdx(~emptyNearestParent)
    obj.data.clusters(nearestParent{p}) = {[obj.data.clusters{nearestParent{p}} pointsIdx(p)]};
end

% Create new cluster IDs
newClusterIDs = cellfun(@sort,obj.data.clusters,'UniformOutput',false);

% Find clusters whose members changed
obj.data.modelIsOutOfDate = ~cellfun(@isequal,clusterIDs,newClusterIDs);

% Restore point order for the clusters that didn't change
obj.data.clusters(~obj.data.modelIsOutOfDate) = clustersOld(~obj.data.modelIsOutOfDate);

% Build null cluster
obj.data.nullCluster = permute(pointsIdx(emptyNearestParent),[2 1]);

disp('Process: All points have been assigned!');

end

