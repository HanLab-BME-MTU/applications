function [iTipVert,iTipEdge] = findTips(edges,nVerts)
%FINDTIPS returns the indices of the skeletal elements which are tips (deg=1)
%
% [iTipVert,iTipEdge] = findTips(edges,nVerts)
%
%Hunter Elliott
%Sometime between 2009 and 2011...
%


%First, get the degree of each vertex and find tips (deg = 1)
vDegree = graphVertDegree(edges,nVerts);
iTipVert = find(vDegree == 1);
nTip = numel(iTipVert);
iTipEdge = zeros(nTip,1);
%Now get the edge associated with each tip
for j = 1:nTip    
    iTipEdge(j) = find(any(bsxfun(@eq,edges,iTipVert(j)),2));    
end