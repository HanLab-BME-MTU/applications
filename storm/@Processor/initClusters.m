function initClusters(obj)

% Initialize clusters
obj.data.clusters = num2cell((1:obj.data.nPoints)');

% Initialize the cluster color
colorTable = hsv(obj.data.nClusters);
colorTable = [colorTable zeros(obj.data.nClusters,1)];
obj.data.clusterColor = colorTable(randperm(obj.data.nClusters),:);

disp('Process: Clusters initialized!');

end