function dissolveClustersSmallerThan(obj,sizeThreshold)

% Find long models
clusterSize = cellfun(@numel,obj.data.clusters);
isBigCluster = clusterSize>=sizeThreshold;

% Update clusters
shortClusters = obj.data.clusters(~isBigCluster);
obj.data.clusters = obj.data.clusters(isBigCluster); 

% Update model lists
obj.removeModels(~isBigCluster);

% Update null cluster
unclustered = permute(horzcat(shortClusters{:}),[2 1]);
obj.data.nullCluster = [obj.data.nullCluster;unclustered];

disp('Process: Small clusters dissolved!');

end


