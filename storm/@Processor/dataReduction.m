function dataReduction(obj,initialEdgeRadius,reductionRuns)
nPointsInitial = obj.data.nPoints;
for i = 1:reductionRuns
    
    % Even number of nodes
    if mod(size(obj.data.points,1),2)
        obj.data.points = obj.data.points(1:end-1,:);
        obj.data.intensity = zeros(obj.data.nPoints,1);
        obj.data.frame = zeros(obj.data.nPoints,1);
    end
    
    % Compute the distance based weights
    edgeLength = obj.initEdges(initialEdgeRadius);
    sigma = initialEdgeRadius/3;
    obj.data.weights = exp(-edgeLength.^2/(2*sigma^2));
    obj.data.weights = reshape(obj.data.weights,length(obj.data.weights),1);
    
    % Find the max weighted matching
    matching = maxWeightedMatching(obj.data.nPoints,obj.data.edges,obj.data.weights);
    
    obj.data.edges = obj.data.edges(matching,:);
    
    fprintf('Process: %d/%d points dropped!\n',obj.data.nPoints-2*nnz(matching),obj.data.nPoints);
    
    % Compute the centers of the point pairs
    obj.data.points = (obj.data.points(obj.data.edges(:,1),:) + obj.data.points(obj.data.edges(:,2),:))/2;
end

obj.data.intensity = zeros(obj.data.nPoints,1);
obj.data.frame = zeros(obj.data.nPoints,1);

fprintf('Process: Data set reduced: %d -> %d\n',nPointsInitial,obj.data.nPoints);
end