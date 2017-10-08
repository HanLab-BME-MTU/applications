function [ S, T, nlms_mip3 ] = getSkeletonFromNLMS( I )
%ANALYZENLMS Produce skeleton using NLMS

gcp;

I = double(I);
F = OrientationSpaceRidgeFilter(1./2/pi./2,[],8,'none');
R = I*F;
R3 = R.getResponseAtOrderFT(3);
[original.maxima,~,original.maximaV] = R.getRidgeOrientationLocalMaxima;
nlms = R3.nonLocalMaximaSuppressionPrecise(original.maxima);


nms = nlms(:,:,1);
nlms_mip3 = nanmax(nlms,[],3);

T = thresholdOtsu(nlms_mip3);
% T = 0;
nlms_binary = nlms_mip3 > T;
nms_binary = nms > T;
% non_nms_binary = nlms_binary & ~nms_binary;


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

