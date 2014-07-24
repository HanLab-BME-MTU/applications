function edgeDir = getEdgePathDir(vertices,edges,edgePaths)
%
%Since I'm an idiot and the direction of the edge-paths is not consistent
%relative to the ordering of the vertices in the edge matrix, this function
%determines which way they go. 
%
% If the edgePath goes in the same direction as the edges entry, there will
% be a 1 in the output vector, if it goes the opposite way it will be a -1
%
% Hunter Elliott
% 2/2012


nEdges = numel(edgePaths);
edgeDir = nan(nEdges,1);
%Go through each edge and check which vertex the starting edgePath point is
%closest to...
for j = 1:nEdges
    
    %Dist from first vertex in edge matrix to first point on edgePath
    dToStart = sqrt(sum((vertices(edges(j,1),:) - edgePaths{j}(1,:)) .^2));
    %Dist from first vertex in edge matrix to last point on edgePath
    dToEnd = sqrt(sum((vertices(edges(j,1),:) - edgePaths{j}(end,:)) .^2));
    
    if dToStart <= dToEnd
        edgeDir(j) = 1;
    else
        edgeDir(j) = -1;
    end    
    
end






