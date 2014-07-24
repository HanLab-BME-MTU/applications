function branchRadii = branchRadii(skelGraph,maskProp)

%Shitty function for getting branch radii... Uses the smoothed surface to
%avoid loading and scaling the mask, calcualting the distance transform and
%sampling it. INstead we use the magic of KDTree.

nEdges = size(skelGraph.edges,1);

branchRadii = cell(nEdges,1);

% %Just the tip, ;) for now....
% [iTipVert,iTipEdge] = findTips(skelGraph.edges,size(skelGraph.vertices,1));
% nTip = numel(iTipVert);
% 

nPtsPerEdge = cellfun(@(x)(size(x,1)),skelGraph.edgePaths);

%Combine all the edge points so we create the KDTree once 
allEdgePts = vertcat(skelGraph.edgePaths{:});

%Now find the closest mask surface point for each edge point
[~,closestDistAll] = KDTreeClosestPoint(maskProp.SmoothedSurface.vertices(:,[2 1 3]),allEdgePts);

%Now separate out the radii into the individual branches
indArray = arrayfun(@(x)(repmat(x,nPtsPerEdge(x),1)),1:nEdges,'Unif',false);
indArray = vertcat(indArray{:});
for j = 1:nEdges
    branchRadii{j} = closestDistAll(indArray  == j);
end





