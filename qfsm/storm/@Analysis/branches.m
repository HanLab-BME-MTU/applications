function branches(obj,myImarisVisualizer,edgeRadius,minAngleThreshold,maxAngleThreshold,maxDistanceThreshold)

% Initialize the edges
pro = Processor(obj.data);
pro.initEdges(edgeRadius);
pro.updateEdges();

% Clear edges leading to a model-less cluster
isBigCluster = obj.data.modelLength > 0;
clustersWithModels = find(isBigCluster);
isModelEdge = ismember(obj.data.edges,clustersWithModels);
isModelEdge = isModelEdge(:,1) & isModelEdge(:,2);
edges = obj.data.edges(isModelEdge,:);
nEdges = size(edges,1);

branchingAngles = [];

% Workaround
branchingPoints = zeros(1,3);

for e=1:nEdges
    clusterIdx1 = edges(e,1);
    clusterIdx2 = edges(e,2);
    
    % Bezier curve end points
    p1s = obj.data.modelBezCP{clusterIdx1}(1,:);
    p1e = obj.data.modelBezCP{clusterIdx1}(end,:);
    p2s = obj.data.modelBezCP{clusterIdx2}(1,:);
    p2e = obj.data.modelBezCP{clusterIdx2}(end,:);
    
    % Distance between the end points an the second curve
    [dist1sOn2,t1sOn2] = distancePointBezier(obj.data.modelBezCP{clusterIdx2},p1s);
    [dist1eOn2,t1eOn2] = distancePointBezier(obj.data.modelBezCP{clusterIdx2},p1e);
    [dist2sOn1,t2sOn1] = distancePointBezier(obj.data.modelBezCP{clusterIdx1},p2s);
    [dist2eOn1,t2eOn1] = distancePointBezier(obj.data.modelBezCP{clusterIdx1},p2e);
    
    % The corresponding projection point on the second curve
    proj1sOn2 = renderBezier(obj.data.modelBezCP{clusterIdx2},t1sOn2);
    proj1eOn2 = renderBezier(obj.data.modelBezCP{clusterIdx2},t1eOn2);
    proj2sOn1 = renderBezier(obj.data.modelBezCP{clusterIdx1},t2sOn1);
    proj2eOn1 = renderBezier(obj.data.modelBezCP{clusterIdx1},t2eOn1);
    
    potentialBranchingPoints = [proj1sOn2+p1s;proj1eOn2+p1e;proj2sOn1+p2s;proj2eOn1+p2e]/2;
    
    dist = [dist1sOn2;dist1eOn2;dist2sOn1;dist2eOn1];
    
    % Normalized tangent vector at the position of the projection
    [~,or1sOn2] = tangentBezier(obj.data.modelBezCP{clusterIdx2},t1sOn2);
    [~,or1eOn2] = tangentBezier(obj.data.modelBezCP{clusterIdx2},t1eOn2);
    [~,or2sOn1] = tangentBezier(obj.data.modelBezCP{clusterIdx1},t2sOn1);
    [~,or2eOn1] = tangentBezier(obj.data.modelBezCP{clusterIdx1},t2eOn1);
    
    % Normalized tangent vector at the end points of the curves
    [~,or1s] = tangentBezier(obj.data.modelBezCP{clusterIdx1},0);
    [~,or1e] = tangentBezier(obj.data.modelBezCP{clusterIdx1},1);
    [~,or2s] = tangentBezier(obj.data.modelBezCP{clusterIdx2},0);
    [~,or2e] = tangentBezier(obj.data.modelBezCP{clusterIdx2},1);
    
    orEndPoints = [or1s;or1e;or2s;or2e];
    orProjections = [or1sOn2;or1eOn2;or2sOn1;or2eOn1];
    
    % Compute branching angles
    angle = acosd(dot(orEndPoints,orProjections,2));
    angle(angle > 90) = 180 - angle(angle > 90);
    
    % Check branching conditions
    isCloseEnough = dist < maxDistanceThreshold;
    hasTheRightAngle = angle < maxAngleThreshold & angle > minAngleThreshold;
    potentialBranchingPointsIdx = (1:4);
    if nnz(isCloseEnough)
        branchingAngles = [branchingAngles;angle(isCloseEnough)];
        if nnz(hasTheRightAngle)
            newBranchingPoints = potentialBranchingPoints(potentialBranchingPointsIdx(isCloseEnough & hasTheRightAngle),:);
            branchingPoints = [branchingPoints;newBranchingPoints];
        end
    end
end

% Display branching angles histogram
hist(branchingAngles,20);
xlabel('Branching angle [degrees]');
ylabel('Number of branching points');

% Workaround
branchingPoints = branchingPoints(2:end,:);

% Display the detected branching points with Imaris
myImarisVisualizer.displayPoints(branchingPoints,20,[1.0 0.0 0.0 0.0]);

disp('Branching: Done!');

end