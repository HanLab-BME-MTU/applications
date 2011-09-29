function [tipPathLengths,tipPaths] = analyzeSkeletonTipPaths(vertices,edges,edgePaths,edgeLabels)




nVert = size(vertices,1);


%Find tips, as these are our starting points
[iTipVerts,iTipEdges] = findTips(edges,nVert);
nTips = numel(iTipVerts);


%Find cell-body points, as these are our end-points. We already have the
%edges, so we need to find the vertices
bodyVerts = find(findBodyVerts(vertices,edges,edgeLabels));
nBodyVerts = numel(bodyVerts);

nEdges = numel(edgePaths);

%Get the length of each edge path.
edgeLen = cellfun(@(x)(sum(sqrt(sum(diff(x,1,1) .^2,2)))),edgePaths);

%And get the weighted adjacency matrix, with the edge lengths as weights
adjMat = sparse(vertcat(edges(:,1),edges(:,2)),vertcat(edges(:,2),edges(:,1)),repmat(edgeLen,[2 1]),nVert,nVert,2*nEdges);


tipPathLengths = nan(nVert,1);
tipPaths = cell(nVert,1);

for j = 1:nTips
    
    for k = 1:nBodyVerts
        
        [currPL(k) currPath{k}] = graphshortestpath(adjMat,iTipVerts(j),bodyVerts(k),'Directed',false);        
        
    end            
    %Keep only the shortest of the routes to the body vertices
    [tipPathLengths(iTipVerts(j)),iShortest] = min(currPL);
    tipPaths{iTipVerts(j)} = currPath{iShortest};
    
end

