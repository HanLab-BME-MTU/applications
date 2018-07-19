function updateEdges(obj)
obj.data.edges = obj.data.parents(obj.data.initialEdges);

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