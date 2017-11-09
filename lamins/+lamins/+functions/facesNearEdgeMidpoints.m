function EF = facesNearEdgeMidpoints(edges,faces)
% function EF = facesNearEdgeMidpoints(edges,faces)
% facesNearEdgeMidpoints finds faces proximal to the middle of each edge
% edges and faces are both connected component structures as output by bwconncomp
    import lamins.functions.*;
    import connectedComponents.*;
    faceL = labelmatrix(faces);
%     faceL = imdilateSafe(faceL,strel('square',3));
    midpoints = edges;
    midpoints.PixelIdxList = num2cell(getEdgeMidpoints(edges));
    midpoints = ccDilate(midpoints,strel('square',5));
    EF = cellfun(@(idx) nonzeros(unique(faceL(idx))),midpoints.PixelIdxList,'UniformOutput',false);
end
