function [ midpoints ] = getEdgeMidpoints( edges )
%getEdgeMidpoints Get the midpoints of edges in a connected compnonents
%structure with a sorted PixelIdxList
    midpoints = cellfun(@(x) single( x( ceil(end/2) ) ),edges.PixelIdxList);
end

