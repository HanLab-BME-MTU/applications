function dataReduction(obj,initialEdgeRadius,reductionRuns)
% function dataReduction(obj,initialEdgeRadius,reductionRuns)
% SYNOPSIS:
% Reduces the total number of points by building pairs of points and
% replacing them with their average. A distance based criterion is used to 
% group the points and the optimal subset is found with a maximum weight 
% matching. Each iteration reduces the number of points by a little bit
% more than a factor 2. Points not part of the matching are discarded.
%
% REQUIRED INPUTS:         
% - initialEdgeRadius
% The maximum distance between two points so that they can still be merged.
% 
% - reductionRuns
% Number of iterations. At each iteration the number of points is reduced
% by a little bit more than a factor 2.
% 
% OPTIONAL INPUTS:
%
% NEEDED PROPERTIES: 
% - obj.data.points
% - obj.data.nPoints
% - obj.data.edges
%
% MODIFIED PROPERTIES:
% - obj.data.points 
% - obj.data.intensity
% - obj.data.frame
% - obj.data.weights
% - obj.data.edges
%
% OUTPUTS:
%
% Pascal Bérard, October 2011

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