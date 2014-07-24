function stop = mergeClustersBIC(obj,maxDegreeBezier,maxCurvature,fitMethod,betaVar,modeVar)

% Compute the edge weights
disp('Process: Computing association weights ...');
edgePointOrder = obj.updateBICWeights(maxDegreeBezier,maxCurvature,fitMethod,betaVar,modeVar);

if numel(obj.data.weights(obj.data.weights>0))==0
    stop = 1;
else
    stop = 0;
        
    % Find the max weighted matching
    matching = maxWeightedMatching(obj.data.nClusters,obj.data.edges,obj.data.weights);
    % matching = maxGreedyMatching(obj.data.edges,obj.data.weights);
    matching(obj.data.weights<0) = false;
    fprintf('Process: Number of matched edges: %u\n',nnz(matching));
    
    obj.data.edges = obj.data.edges(matching,:);
    obj.data.weights = obj.data.weights(matching,:);
    edgePointOrder = edgePointOrder(matching,:);
    
    % Merge the clusters
    mergedClusters = arrayfun(@(a,b) [obj.data.clusters{a} obj.data.clusters{b}],obj.data.edges(:,1),obj.data.edges(:,2),'UniformOutput',false);   
    mergedClusters = cellfun(@(a,b) a(b),mergedClusters,edgePointOrder,'UniformOutput',false);
    
    mergedClustersColor = obj.data.clusterColor(obj.data.edges(:,1),:);
    
    % Add unmatched clusters to cluster list
    nClusters = size(obj.data.clusters,1);
    unmatchedClustersIdxs = setdiff((1:nClusters)',[obj.data.edges(:,1); obj.data.edges(:,2)]);
    unmatchedClustersColor = obj.data.clusterColor(unmatchedClustersIdxs,:);
        
    obj.data.clusters = [mergedClusters;obj.data.clusters(unmatchedClustersIdxs)];    
    obj.data.clusterColor = [mergedClustersColor;unmatchedClustersColor];

    obj.data.modelType = [obj.edgeModelType(matching);obj.data.modelType(unmatchedClustersIdxs)];
    obj.data.modelLength = [obj.edgeModelLength(matching);obj.data.modelLength(unmatchedClustersIdxs)];
    obj.data.modelRes = [obj.edgeModelRes(matching);obj.data.modelRes(unmatchedClustersIdxs)];
    obj.data.modelProj = [obj.edgeModelProj(matching);obj.data.modelProj(unmatchedClustersIdxs)];
    obj.data.modelBezCP = [obj.edgeModelBezCP(matching);obj.data.modelBezCP(unmatchedClustersIdxs)];
    obj.data.modelVar = [obj.edgeModelVar(matching);obj.data.modelVar(unmatchedClustersIdxs)];
    
    obj.data.modelIsOutOfDate = [false(nnz(matching),1);obj.data.modelIsOutOfDate(unmatchedClustersIdxs)];    
    
    disp('Process: Cluster merged!');  
end
end
