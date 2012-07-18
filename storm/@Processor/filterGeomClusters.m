function filterGeomClusters(obj,filterData,nSigmaDist)

% Reference data set non-noise points
refPointIdx = horzcat(filterData.clusters{:});
refPoints = filterData.points(refPointIdx,:);
refParents = filterData.parents(refPointIdx);
refError = filterData.error(refPointIdx,:);
refScaPoints = refPoints ./ refError;
refModelBezCP = filterData.modelBezCP;
refModelType = filterData.modelType;

% Normalize points
scaPoints = obj.data.points ./ obj.data.error; % Scaled points

% Find the neighbor points in the reference data set
nei = KDTreeBallQuery(refScaPoints,scaPoints,100); % Neighbors

% Find the corresponding parent models
neiPar = cellfun(@(a) refParents(a(:)),nei,'UniformOutput',false); % Neighbor parent model
neiPar = cellfun(@unique,neiPar,'UniformOutput',false); % Remove double entries

% Compute the distance to these models
obj_data_error = obj.data.error;
filterPoint = true(obj.data.nPoints,1);

parfor i=1:obj.data.nPoints
    dist = arrayfun(@(b) distancePointBezier(refModelBezCP{b}./repmat(obj_data_error(i,:),refModelType(b)+1,1),scaPoints(i,:)),neiPar{i});
        
    if min(dist) < nSigmaDist
        filterPoint(i) = false;
    end
end

% If a cluster contains one of these points then dissolve the cluster
filterPointIdx = find(filterPoint);
clusterToDissolveIdx = cellfun(@(a) all(ismember(a,filterPointIdx)),obj.data.clusters);

% Update clusters
clusters = obj.data.clusters(~clusterToDissolveIdx);
clustersToDissolve = obj.data.clusters(clusterToDissolveIdx);
unclustered = horzcat(clustersToDissolve{:})';
obj.data.clusters = clusters;
obj.data.nullCluster = [obj.data.nullCluster;unclustered];

end




