function mergeClustersGeom(obj,modelLength,angleThreshold)

% Compute the meanLocPrec
meanLocPrec = mean(obj.data.error,1);

% Compute the edge weights
[obj.data.weights,hausDorffEdges] = obj.updateGeomWeights(meanLocPrec,modelLength,angleThreshold);

% Remove disabled edges
obj.data.edges = obj.data.edges(obj.data.weights>=0,:);
hausDorffEdges = hausDorffEdges(obj.data.weights>=0,:);
obj.data.weights = obj.data.weights(obj.data.weights>=0);

% Find the max weighted matching
matching = maxWeightedMatching(size(obj.data.clusters,1),obj.data.edges,obj.data.weights);
obj.data.edges = obj.data.edges(matching,:);
hausDorffEdges = hausDorffEdges(matching,:);
obj.data.weights = obj.data.weights(matching,:);
disp(['Process: Number of matched edges: ' num2str(nnz(matching))]);

% Merge the clusters
mergedClusters = arrayfun(@(a,b) [obj.data.clusters{a} obj.data.clusters{b}],obj.data.edges(:,1),obj.data.edges(:,2),'UniformOutput',false);
mergedClustersColor = obj.data.clusterColor(obj.data.edges(:,1),:);

% Add unmatched clusters to cluster list
unmatchedClustersIdxs = setdiff((1:obj.data.nClusters)',[obj.data.edges(:,1); obj.data.edges(:,2)]);
unmatchedClustersColor = obj.data.clusterColor(unmatchedClustersIdxs,:);

obj.data.clusters = [mergedClusters;obj.data.clusters(unmatchedClustersIdxs)];
obj.data.clusterColor = [mergedClustersColor;unmatchedClustersColor];

% Replace the edges with the indices of the hausDorffEdges so that they can
% be stored with the saveEdgesToHistory method and displayed with
% edgesHistoryLast.
obj.data.edges = hausDorffEdges;

disp('Process: Cluster merged!');

end