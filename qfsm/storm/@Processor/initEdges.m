function edgeLength = initEdges(obj,initialEdgeRadius)

% Find the neighbors of all the points
radii = repmat(initialEdgeRadius,obj.data.nPoints,1);
[neighbors,distances] = KDTreeBallQuery(obj.data.points,obj.data.points,radii);

% Store the result of the kd-tree query
obj.data.neighbors = cellfun(@(a) a(2:end),neighbors,'UniformOutput',false);

% Transform the cellarray into a vector
nNeighbors = cellfun(@numel,neighbors);
lastNeighbor = cumsum(nNeighbors);
centers = (1:size(lastNeighbor,1))';
neighborsVector = vertcat(neighbors{:});

% Construct the vector with the corresponding centers
centers = arrayfun(@(a,b) ones(b,1)*a,centers,nNeighbors,'UniformOutput',false);
centers = vertcat(centers{:});

% Create edge list
edges = [neighborsVector centers];

% Put the smaller cluster index to the left [1 2] and [2 1] => [1 2]
edges = sort(edges,2);

% Remove the edge duplicates
[edges,uniqueIdx] = unique(edges,'rows');

% Remove self-edges
nonSelfEdgeIdx = edges(:,1)~=edges(:,2);
edges = edges(nonSelfEdgeIdx,:);

% Return the length of the edges
if nargout > 0
    distancesVector = vertcat(distances{:});
    distancesVector = distancesVector(uniqueIdx);
    edgeLength = distancesVector(nonSelfEdgeIdx);
end

% Initialize the weights
obj.data.weights = -ones(size(edges,1),1);

% Save as initial edges
obj.data.initialEdges = edges;
obj.data.edges = edges;

disp('Process: Initial edges created!');

end
