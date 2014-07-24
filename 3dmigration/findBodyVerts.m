function vertLabels = findBodyVerts(vertices,edges,edgeLabels)



nVerts = size(vertices,1);

vertLabels = zeros(nVerts,1);

for j = 1:nVerts
    
    isEdge = any(edges == j,2);
    
    vertLabels(j) = nnz(edgeLabels(isEdge) == 2);
    
    
end




