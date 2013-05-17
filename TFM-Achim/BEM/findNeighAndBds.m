function [neigh,bounds,bdPtsID]=findNeighAndBds(p,t)
% generate neighbors, bounding box for each basis function and the convex
% hull of the mesh.
% This function has been tested for square and hexagonal lattices.
numPts=size(p,1);
neigh(numPts).cand=[];
bounds(numPts).x=[];
bounds(numPts).y=[];

for idx=1:length(p(:,1))
    [row,~] = find(t==idx);
    cand=setdiff(unique(t(row,:)),idx);
    neigh(idx).cand=round(cand(:)');
    neigh(idx).pos =p(neigh(idx).cand,:);
    bounds(idx).x  =[min(neigh(idx).pos(:,1)) max(neigh(idx).pos(:,1))];
    bounds(idx).y  =[min(neigh(idx).pos(:,2)) max(neigh(idx).pos(:,2))];
end
bdPtsID=convhull(p(:,1),p(:,2));