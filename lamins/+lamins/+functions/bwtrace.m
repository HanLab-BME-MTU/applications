function [ edges_cc, vertices_cc, pt_edge_pairs ] = bwtrace( bw, pts )
%bwtrace 
%%
% Also see skel2graph

pts_label = uint16(pts);
Npts = sum(pts_label(:));
pts_label(pts) = 1:Npts;
dilated_pts = imdilate(pts_label,strel('square',3));

edges = bw & ~dilated_pts;
edges_cc = bwconncomp(edges);
edges_lm = labelmatrix(edges_cc);
edges_lm = uint16(imdilate(edges_lm,strel('square',3))).*uint16(bw);
rp = regionprops(edges_lm,'PixelIdxList');
edges_cc.PixelIdxList = {rp.PixelIdxList};
edges_cc.lengths = cellfun(@length,edges_cc.PixelIdxList);

nonzero = dilated_pts & edges_lm;
pt_edge_pairs = unique([dilated_pts(nonzero) edges_lm(nonzero)],'rows');

for i=1:edges_cc.NumObjects
    edges_cc.vertices{i} = pt_edge_pairs( pt_edge_pairs(:,2) == i, 1);
end

vertices_cc.Connectivity = 0;
vertices_cc.ImageSize = edges_cc.ImageSize;
vertices_cc.NumObjects = Npts;
vertices_cc.PixelIdxList = num2cell(find(pts_label))';

for i=1:Npts
    vertices_cc.edges{i} = pt_edge_pairs( pt_edge_pairs(:,1) == i, 2);
end


end

