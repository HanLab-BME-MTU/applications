function [S, TfacesMapDilatedMin] = rethresholdSkeletonFromNLMS(S,I,nlms_mip3,nms_binary,T)
%rethresholdSkeletonFromNLMS Use faces from Skeleton to recalculate
%regional thresholds
    if(isscalar(T))
        T = T*ones(size(nlms_mip3));
    end
    Tfaces = zeros(S.faces.NumObjects,1);
    TfacesMap = zeros(size(nlms_mip3));
    for f=1:S.faces.NumObjects
        try
            Tfaces(f) = thresholdOtsu(nlms_mip3(S.faces.PixelIdxList{f}));
            TfacesMap(S.faces.PixelIdxList{f}) = Tfaces(f);
        catch
        end;
    end;

    TfacesMapDilated = imdilate(TfacesMap,strel('disk',3));

    TfacesMapDilated(TfacesMapDilated == 0) = T(TfacesMapDilated == 0);
    TfacesMapDilatedMin = min(TfacesMapDilated,T);

    nlms_binary = nlms_mip3 > TfacesMapDilatedMin;

    % From getSkeletonFromNLMS
    nms_skel = bwmorph(nms_binary,'skel',Inf);
    nms_skel_thresh = nms_skel.* nlms_binary;
    nms_skel_thresh_wo_bp = lamins.functions.bwRemoveBranchPoints(nms_skel_thresh);
    nlms_fragments = nlms_binary & ~nms_skel_thresh_wo_bp;
    nlms_fragments_cc = bwconncomp(nlms_fragments);
    nms_skel_wo_bp_cc = bwconncomp(nms_skel_thresh_wo_bp);
    bridge = lamins.functions.minimalBridge(nlms_fragments_cc,nms_skel_wo_bp_cc,I);
    nms_skel_bridged = nms_skel | bridge;
    nms_skel_bridged_skel = bwmorph(nms_skel_bridged,'skel',Inf);
    S = lamins.classes.Skeleton(nms_skel_bridged_skel);
    %  lamins.classes.Skeleton.cleanup without removing edges
    S.convertShortEdgesToVertices(2);
    S.reduceVerticesToPoints;
    FE = S.faceEdges;
    f = find(cellfun(@length,FE) == 0);
    S.deleteFaces(f);
end