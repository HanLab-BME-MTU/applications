function dissolveModelsLessDenseThan(obj,densityThreshold)

% Find dense clusters
clusterSize = cellfun(@numel,obj.data.clusters);
modelLength = obj.data.modelLength;
modelLength(modelLength == 0) = -1;
isDenseCluster = (clusterSize./modelLength)>=densityThreshold;

% Update clusters
shortClusters = obj.data.clusters(~isDenseCluster);
obj.data.clusters = obj.data.clusters(isDenseCluster);

% Update model lists
obj.removeModels(~isDenseCluster);

% Update null cluster
unclustered = permute(horzcat(shortClusters{:}),[2 1]);
obj.data.nullCluster = [obj.data.nullCluster;unclustered];

disp('Process: Sparse models dissolved!');

end


