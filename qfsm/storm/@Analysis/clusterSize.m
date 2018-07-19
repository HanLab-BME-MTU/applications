function clusterSize(obj)
nBins = round(sqrt(obj.data.nClusters));
nPerModel = cellfun(@numel,obj.data.clusters);

% Display the cluster size histogram
hist(nPerModel,nBins);
xlabel('Cluster size');
ylabel('Number of clusters');
end