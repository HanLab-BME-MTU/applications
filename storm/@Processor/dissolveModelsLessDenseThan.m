function dissolveModelsLessDenseThan(obj,densityThreshold)

% Find dense clusters
isDenseCluster = (obj.data.clusterSize./obj.data.modelLength)>=densityThreshold;

% Update clusters
shortClusters = obj.data.clusters(~isDenseCluster);
obj.data.clusters = obj.data.clusters(isDenseCluster);

% Update model lists
obj.removeModels(~isDenseCluster);

% Update null cluster
unclustered = horzcat(shortClusters{:})';
obj.data.nullCluster = [obj.data.nullCluster;unclustered];

disp('Process: Sparse models dissolved!');

end


