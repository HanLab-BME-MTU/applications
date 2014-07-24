function bundles(obj,imaris,edgeRadius,distanceThreshold,minOverlapThreshold)

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

isBundle = false(nEdges,2);
for e=1:nEdges
    clusterIdx1 = edges(e,1);
    clusterIdx2 = edges(e,2);
    
    % Model end points
    p1s = obj.data.modelBezCP{clusterIdx1}(1,:);
    p1e = obj.data.modelBezCP{clusterIdx1}(end,:);
    p2s = obj.data.modelBezCP{clusterIdx2}(1,:);
    p2e = obj.data.modelBezCP{clusterIdx2}(end,:);
    
    % Bezier curve control points
    cP1 = obj.data.modelBezCP{clusterIdx1};
    cP2 = obj.data.modelBezCP{clusterIdx2};
    
    % Distance of the endpoints to the other curve
    [d1s,t1s] = distancePointBezier(cP2,p1s);
    [d1e,t1e] = distancePointBezier(cP2,p1e);
    [d2s,t2s] = distancePointBezier(cP1,p2s);
    [d2e,t2e] = distancePointBezier(cP1,p2e);
    
    % Check bundle conditions for model 2
    if abs(t1s-t1e) > minOverlapThreshold
        if t1s == 0
            d1s = distancePointBezier(cP1,cP2(1,:));
        elseif t1s == 1
            d1s = distancePointBezier(cP1,cP2(end,:));
        end
        if t1e == 0
            d1e = distancePointBezier(cP1,cP2(1,:));
        elseif t1e == 1
            d1e = distancePointBezier(cP1,cP2(end,:));
        end
        if d1s < distanceThreshold && d1e < distanceThreshold
            isBundle(e,2) = true;
        end
    end
    
     % Check bundle conditions for model 1
    if abs(t2s-t2e) > minOverlapThreshold
        if t2s == 0
            d2s = distancePointBezier(cP2,cP1(1,:));
        elseif t2s == 1
            d2s = distancePointBezier(cP2,cP1(end,:));
        end
        if t2e == 0
            d2e = distancePointBezier(cP2,cP1(1,:));
        elseif t2e == 1
            d2e = distancePointBezier(cP2,cP1(end,:));
        end
        if d2s < distanceThreshold && d2e < distanceThreshold
            isBundle(e,1) = true;
        end
    end
end

edgesVec = edges(:);
isBundleVec = isBundle(:);
modelsInBundle = edgesVec(isBundleVec);
modelsInBundle = unique(modelsInBundle)
clustersInBundle = obj.data.clusters(modelsInBundle);
pointsInBundle = obj.data.points(horzcat(clustersInBundle{:})',:);

% Display detected bundles in Imaris
if ~isempty(pointsInBundle)
    imaris.displayPoints(pointsInBundle,5,[0.0 0.0 1.0 0.0]);
end

disp('Bundling: Done!');

end
