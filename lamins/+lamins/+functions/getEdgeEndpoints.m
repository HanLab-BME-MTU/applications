function [ endpts ] = getEdgeEndpoints( edges )
%getEdgeEndpts Get the midpoints of edges in a connected compnonents
%structure with a sorted PixelIdxListSummary of this function goes here

    endpts = cellfun(@(x) [ x(1) x(end) ],edges.PixelIdxList,'UniformOutput',false);
    endpts = vertcat(endpts{:});
end

