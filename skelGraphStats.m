function gStats= skelGraphStats(vertices,edges,edgePaths)
%UNDER CONSTRUCTION -HLE




if nargin < 3 || isempty(vertices) || isempty(edges) || isempty(edgePaths)
    error('You must input the full skeleton graph, including vertices, edges and edgePaths!')
end


%MORE INPUT CHECKING TEMP - HLE

nVert = size(vertices,1);
nEdges = size(edges,1);

gStats.vertexDegree = zeros(max(edges(:)),1);

%Get the degree of each vertex
for j = 1:nEdges
    
    %Make sure it's not a spur
    if ~any(edges(j,:)==0)
        %Add this edge to the degree count of each vertex it connects
        gStats.vertexDegree(edges(j,:)) = gStats.vertexDegree(edges(j,:)) + 1;
    end
end

