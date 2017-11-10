function [edges_cc,vertices_cc, pairs, skel] = getSkeletonGraph(obj,varargin)
    skel = obj.getSkeleton(varargin{:});
    branch_pts = bwmorph(skel,'branchpoints');
    end_pts = bwmorph(skel,'endpoints');
    vertices = branch_pts | end_pts;
    [edges_cc, vertices_cc, pairs] = bwtrace(skel,vertices);
end

