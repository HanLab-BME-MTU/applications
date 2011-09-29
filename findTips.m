function [iTipVert,iTipEdge] = findTips(edges,nVerts)

%First, get the degree of each vertex and find tips (deg = 1)
vDegree = graphVertDegree(edges,nVerts);
iTipVert = find(vDegree == 1);
nTip = numel(iTipVert);
iTipEdge = zeros(nTip,1);
%Now get the edge associated with each tip
for j = 1:nTip    
    iTipEdge(j) = find(any(bsxfun(@eq,edges,iTipVert(j)),2));    
end