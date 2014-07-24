function saveClustersToHistory(obj)
obj.data.clustersHistory{numel(obj.data.clustersHistory)+1} = obj.data.clusters;
end