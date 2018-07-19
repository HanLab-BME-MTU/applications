function updateEdgesEndPoints(obj,dMax)

% Get the model end points (first and last control point)
samplePnts = cellfun(@(a) [a(1,:);a(end,:)],obj.data.modelBezCP,'UniformOutput',false);
num = 2*ones(size(samplePnts,1),1);
last = cumsum(num);
first = last-num+1;
samplePntsVector = vertcat(samplePnts{:});

% Find the neighbor points of the end points
neighPntIdxVector = KDTreeBallQuery(obj.data.points,samplePntsVector,repmat(dMax,size(samplePntsVector,1)));
neighPntIdx = arrayfun(@(a,b) vertcat(neighPntIdxVector{a:b}),first,last,'UniformOutput',false);

% Remove duplicate neighbor point indices
neighPntIdx = cellfun(@unique,neighPntIdx,'UniformOutput',false);
neighPntIdxVector = vertcat(neighPntIdx{:});

% Find the parent of the neighbor points
neighModelIdxVector = obj.data.parents(neighPntIdxVector);

% Build the edge list
nNeighbors = cellfun(@numel,neighPntIdx);
modelIdx = arrayfun(@(a,b) ones(a,1)*b,nNeighbors',1:obj.data.nClusters,'UniformOutput',false);
modelIdxVector = vertcat(modelIdx{:});
obj.data.edges = [modelIdxVector neighModelIdxVector];

% Put the smaller cluster index to the left [1 2] and [2 1] => [1 2]
obj.data.edges = sort(obj.data.edges,2);

% Remove the edge duplicates
obj.data.edges = unique(obj.data.edges,'rows');

% Remove self-edges
obj.data.edges = obj.data.edges(obj.data.edges(:,1)~=obj.data.edges(:,2),:);

% Initialize the weights
obj.data.weights = -ones(size(obj.data.edges,1),1);

disp('Process: Edges updated!');
end