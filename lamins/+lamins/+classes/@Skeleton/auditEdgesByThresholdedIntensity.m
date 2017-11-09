function S = auditEdgesByThresholdedIntensity(S,I)
    F = I.flattenIntensity;
%     T = thresholdOtsu(F);
    % use mask threshold for evaluating region properties
    T_mask = I.maskThresh(F);
    % use otsu threshold for evaluting T_edge
    T_otsu = thresholdOtsu(F);
    B_mask = F > T_mask;
    B_otsu = F > T_otsu;
    T_edge = sum(B_otsu(I.mask(:)))./sum(I.mask(:));
    rp = regionprops(S.edges,double(B_mask),'MeanIntensity');
    S.deleteEdges([rp.MeanIntensity] < T_edge);
end