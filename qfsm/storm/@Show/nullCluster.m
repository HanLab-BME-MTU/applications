function nullCluster(obj)

% Display null cluster
if numel(obj.data.nullCluster)
    obj.imaris.displayPoints(obj.data.points(obj.data.nullCluster,:),obj.pointSize,obj.nullClusterColor,'Display: Null Cluster');
end

end



