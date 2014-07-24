function [weights,hausDorffEdges] = updateGeomWeights(obj,meanLocPrec,modelLength,angleThreshold)

%
% Weights based on the segment distance and the angle between the segments
%

% Copy data
obj_data_edges = obj.data.edges;
obj_data_clusters = obj.data.clusters;

% Init variables
weights = zeros(size(obj_data_edges,1),1);
hausDorffEdges = zeros(size(obj_data_edges,1),2);

% Loop through all the outer edges
parfor i=1:size(obj_data_edges,1)
    % Get the cluster indices
    e1 = obj_data_edges(i,1); % First edge node
    e2 = obj_data_edges(i,2);
    
    % Get the point indices
    c1 = obj_data_clusters{e1}; % First cluster
    c2 = obj_data_clusters{e2};
    
    % Array with all the inner edges between c1 and c2
    [c1Grid,c2Grid] = meshgrid(c1,c2);
    
    % Compute all the weights of the inner edges
    innerWeights = arrayfun(@(a,b) obj.geomClusterWeight(a,b,meanLocPrec,modelLength,angleThreshold),c1Grid,c2Grid);
    
    % Find the inner edge corresponding to the maximum weight
    [rowVal,rowIdx] = max(innerWeights,[],1);
    [colVal,colIdx] = max(rowVal);
    weights(i,1) = colVal;
    
    % Save the inner edge that will be visualized
    hausDorffEdges(i,:) = [c1Grid(rowIdx(colIdx),colIdx) c2Grid(rowIdx(colIdx),colIdx)];
end

end