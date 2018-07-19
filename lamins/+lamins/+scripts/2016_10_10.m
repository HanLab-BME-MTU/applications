I = imread('example.tif');
F = OrientationSpaceRidgeFilter(1./2/pi./2,[],8,'none');
R = F*I2;
I2 = I((-100:100)+586+64,(-100:100)+307+64)
R = F*I2;
[maxima] = R.getRidgeOrientationLocalMaxima;
R3 = R.getResponseAtOrderFT(3);
nlms = R3.nonLocalMaximaSuppressionPrecise(maxima)
nms_skel = bwmorph(nlms(:,:,1) > 0,'skel',Inf);
nms_skel_wo_bp = lamins.functions.bwRemoveBranchPoints(nms_skel);
nms_skel_wo_bp_cc = bwconncomp(nms_skel_wo_bp);
nms_skel_wo_bp_label = labelmatrix(nms_skel_wo_bp_cc);
imshow(label2rgb(nms_skel_wo_bp_label,'hsv','k','shuffle')
imshow(label2rgb(nms_skel_wo_bp_label,'hsv','k','shuffle'),[])
imshow(I2,[]);
imshow(nanmax(nlms,[],3),[]);
imshow(nanmax(nms_skel,[],3),[]);
imshow(nanmax(nms_skel_wo_bp,[],3),[]);
imshow(label2rgb(nms_skel_wo_bp_label,'hsv','k','shuffle'),[])
imshow(I2,[]);
imshow(nanmax(nlms,[],3),[]);
imshow(nanmax(nms_skel,[],3),[]);
imshow(nanmax(nms_skel_wo_bp,[],3),[]);
nlms_mip = nanmax(nlms,[],3);
imshow(nlms_mip > 0,[]);
nlms_mip_skel = bwmorph(nmls_mip,'skel',Inf);
nlms_mip_skel = bwmorph(nlms_mip,'skel',Inf);
imshow(nlms_mip,[]);
imshow(nlms_mip > 0,[]);
imshow(nlms_mip,[]);
imshow(nlms_mip > 0,[]);
imshow(nlms_mip,[]);
thresholdOtsu(nlms_mip)
imshow(nlms_mip > 86,[]);
figure; imshow(nms_skel.*nlms(:,:,1),[]);
figure; imshow(nms_skel.*nlms(:,:,1) > 86,[]);
nms_skel_thresh = nms_skel.*nlms(:,:,1) > 86;
nms_mip_thresh = nlms_mip > 86;
nlms_mip_thresh = nlms_mip > 86;
imshowpair(nms_mip_thresh,nlms_mip_thresh)
imshowpair(nms_skel_thresh,nlms_mip_thresh)
figure; imshow(nms_skel.*nlms(:,:,1) > 86,[]);
nms_skel_thresh
nms_skel_thresh_wo_bp = lamins.functions.bwRemoveBranchPoints(nms_skel_thresh);
nms_skel_thresh_wo_bp_cc = bwconncomp(nms_skel_thresh_wo_bp);
nms_skel_thresh_wo_bp_label = labelmatrix(nms_skel_thresh_wo_bp_cc);
imshowpair(nms_skel_thresh,nms_skel_thresh_wo_bp)
imshow(label2rgb(nms_skel__thresh_wo_bp_label,'hsv','k','shuffle'),[])
imshow(label2rgb(nms_skel_thresh_wo_bp_label,'hsv','k','shuffle'),[])
imshowpair(nms_skel_thresh_wo_bp,nlms_mip_thresh)
nlms_mip_thresh & ~nms_skel_thresh_wo_bp
imshow(nlms_mip_thresh & ~nms_skel_thresh_wo_bp,[]);
nlms_fragments = nlms_mip_thresh & ~nms_skel_thresh_wo_bp;
nlms_fragments_cc = bwconncomp(nlms_fragments);
imshow(label2rgb(labelmatrix(nlms_fragments_cc),'hsv','k','shuffle'),[])
imshow(label2rgb(labelmatrix(connectedComponents.ccDilate(nlms_fragments_cc)),'hsv','k','shuffle'),[])
imshow(label2rgb(labelmatrix(connectedComponents.ccDilate(nlms_fragments_cc,ones(3))),'hsv','k','shuffle'),[])
connectedComponents.ccDilate(nlms_fragments_cc,ones(3))
nlms_fragments_cc2 = connectedComponents.ccDilate(nlms_fragments_cc,ones(3));
cellfun(@(x) unique(nms_skel_thresh_wo_bp_label(x)),nlms_fragments_cc2)
cellfun(@(x) unique(nms_skel_thresh_wo_bp_label(x)),nlms_fragments_cc2.PixelIdxList)
cellfun(@(x) unique(nms_skel_thresh_wo_bp_label(x)),nlms_fragments_cc2.PixelIdxList,'Unif',false)
fragment_segments = cellfun(@(x) unique(nms_skel_thresh_wo_bp_label(x)),nlms_fragments_cc2.PixelIdxList,'Unif',false);
fragment_segments = cellfun(@(x) x(x~=0),fragment_segments);
fragment_segments = cellfun(@(x) x(x~=0),fragment_segments,'Unif',false);
fragment_segments
doc cellfun
cellfun('length',fragment_segments)
cellfun('length',fragment_segments) > 1
fragment_segments_filt = fragment_segments(cellfun('length',fragment_segments) > 1);
fragment_segments_filt
cc3 = cc2;
nlms_fragments_cc3 = nlms_fragments_cc2;
nlms_fragments_cc3
nlms_fragments_cc3.PixelIdxList = fragment_segments_filt;
nlms_fragments_cc3
nlms_fragments_cc3.PixelIdxList = fragment_segments_filt';
nlms_fragments_cc3
nlms_fragments_cc3.NumObjects = 303;
imshow(nlms_fragments,[]);
nlms_fragments_cc3
labelmatrix(nlms_fragments_cc3)
imshow(labelmatrix(nlms_fragments_cc3),[])
nlms_fragments_cc3.PixelIdxList
nlms_fragments
imshow(nlms_fragments,[])
nlms_fragments_cc3 = nlms_fragments_cc2;
nlms_fragments
nlms_fragments_cc
nlms_fragment_segments
nlms_fragments_filt = nlms_fragments_cc(cellfun('length',fragment_segments) > 1);
nlms_fragments_cc
fragment_segments
nlms_fragments_cc
nlms_fragments_filt = nlms_fragments_cc.PixelIdxList(cellfun('length',fragment_segments) > 1);
nlms_fragments_filt
nlms_fragments_cc3.PixelIdxList = nlms_fragments_filt;
nlms_fragments_cc3
nlms_fragments_cc3.NumObjects = 303;
imshow(labelmatrix(nlms_fragments_cc3,[]);
imshow(labelmatrix(nlms_fragments_cc3,[]));
imshow(labelmatrix(nlms_fragments_cc3),[]);
imshow(label2rgb(labelmatrix(nlms_fragments_cc3)),[]);
imshowpair(labelmatrix(nlms_fragments_cc3) > 0,nms_skel_thresh_wo_bp);
connectedComponents.ccDilate(nms_skel_thresh_wo_bp_cc)
connectedComponents.ccDilate(nms_skel_thresh_wo_bp_cc,ones(3))
nms_skel_thresh_wo_bp_cc_dilated = connectedComponents.ccDilate(nms_skel_thresh_wo_bp_cc,ones(3));
nms_skel_thresh_wo_bp_cc_dilated
nms_skel_thresh_wo_bp_cc_dilated.PixelIdxList
for i=1:439; nms_skel_thresh_wo_bp_cc_dilated.PixelIdxList{i}; end;
count = zeros(201);
for i=1:439; count(nms_skel_thresh_wo_bp_cc_dilated.PixelIdxList{i}) = count(nms_skel_thresh_wo_bp_cc_dilated.PixelIdxList{i})+1 ; end;
imshow(count,[]);
imshow(label2rgb(labelmatrix(nms_skel_thresh_wo_bp_cc_dilated)),[]);
imshow(count,[]);
nlms_fragments
nlms_fragments_cc
regionprops(nlms_fragments_cc,count,'Max')
maxCount = regionprops(nlms_fragments_cc,count,'Max');
maxCount
maxCount.MaxIntensity
[maxCount.MaxIntensity]
[maxCount.MaxIntensity] > 1
nlms_fragments([maxCount.MaxIntensity] > 1)
nlms_fragments_cc.PixelIdxList([maxCount.MaxIntensity] > 1)
[maxCount.MaxIntensity]
maxCount = regionprops(nlms_fragments_cc3,count,'Max');
maxCount
maxCount.MaxIntensity
save('2016_10_10.mat')