function mergeWith(obj,dataToMerge)

offset = dataToMerge.nPoints;

% ROI
obj.roiPosition = [];
obj.roiSize = [];

% Points
obj.rawPoints = [obj.rawPoints; dataToMerge.rawPoints];
obj.intensity = [obj.intensity; dataToMerge.intensity];
obj.frame = [obj.frame; dataToMerge.frame];
obj.error = [obj.error; dataToMerge.error];
obj.points = [obj.points; dataToMerge.points];

obj.orientation = [obj.orientation; dataToMerge.orientation];
obj.magnitude = [obj.magnitude; dataToMerge.magnitude];

neighborsToMerge = cellfun(@(a) a+offset, dataToMerge.neighbors, 'UniformOutput', false);
obj.neighbors = [obj.neighbors; neighborsToMerge];

% Clusters
clustersToMerge = cellfun(@(a) a+offset, dataToMerge.clusters, 'UniformOutput', false);
obj.clusters = [obj.clusters; clustersToMerge];

nullClusterToMerge = dataToMerge.nullCluster+offset;
obj.nullCluster = [obj.nullCluster; nullClusterToMerge];

obj.clusterColor = [obj.clusterColor; dataToMerge.clusterColor];

% Cluster models
obj.modelIsOutOfDate = [obj.modelIsOutOfDate; dataToMerge.modelIsOutOfDate];
obj.modelType = [obj.modelType; dataToMerge.modelType];
obj.modelLength = [obj.modelLength; dataToMerge.modelLength];
obj.modelRes = [obj.modelRes; dataToMerge.modelRes];
obj.modelProj = [obj.modelProj; dataToMerge.modelProj];
obj.modelBezCP = [obj.modelBezCP; dataToMerge.modelBezCP];
obj.modelVar = [obj.modelVar; dataToMerge.modelVar];

% Edges
initialEdgesToMerge = dataToMerge.initialEdges+offset;
obj.initialEdges = [obj.initialEdges; initialEdgesToMerge];

edgesToMerge = dataToMerge.edges+offset;
obj.edges = [obj.edges; edgesToMerge];

obj.weights = [obj.weights; dataToMerge.weights];

% Running time
obj.runTime = [];

obj.edgesHistory = [];
obj.clustersHistory = [];
obj.modelsHistory = [];

obj.simModelBezCP = [obj.simModelBezCP; dataToMerge.simModelBezCP];

disp('Process: Data sets merged!');

end

