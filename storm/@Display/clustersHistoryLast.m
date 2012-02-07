function clustersHistoryLast(obj)
nClusters = numel(obj.data.clustersHistory{end});
colorTable = hsv(nClusters);
colorTable = [colorTable zeros(nClusters,1)];
clusterColor = colorTable(randperm(nClusters),:);
for c=1:nClusters
    points = obj.data.points(obj.data.clustersHistory{end}{c},:);
    obj.imaris.displayPoints(points,4,clusterColor(c,:),'Display: Clusters History');
end
end