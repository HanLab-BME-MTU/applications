function clusters(obj,idx)

% Display a single cluster
if nargin > 1
    points = obj.data.points(obj.data.clusters{idx},:);
    obj.imaris.displayPoints(points,obj.pointSize*1.01,obj.data.clusterColor(idx,:),'Display: Cluster');
else
    
    % Display clusters
    for c=1:obj.data.nClusters
        points = obj.data.points(obj.data.clusters{c},:);
        obj.imaris.displayPoints(points,obj.pointSize,obj.data.clusterColor(c,:),'Display: Clusters');
    end
    
    % Display null cluster
    obj.nullCluster();
    
end

end


