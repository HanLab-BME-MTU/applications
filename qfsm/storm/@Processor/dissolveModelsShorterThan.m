function dissolveModelsShorterThan(obj,thresholdLength)

% Find long models
isLongCluster = obj.data.modelLength>=thresholdLength;

% Update clusters
shortClusters = obj.data.clusters(~isLongCluster);
obj.data.clusters = obj.data.clusters(isLongCluster);

% Update model lists
obj.removeModels(~isLongCluster);

% Update null cluster
unclustered = permute(horzcat(shortClusters{:}),[2 1]);
obj.data.nullCluster = [obj.data.nullCluster;unclustered];

disp('Process: Short models dissolved!');

end


