function [mask,cutThreshHL] = refineEdgeWithSteerableFilter(mask0,image,...
    doPlot,closureRadius)
%REFINEEDGEWITHSTEERABLEFILTER refines cell edge using intensity gradients obtained from a steerable line filter
%
%SYNOPSIS [mask,cutThreshHL] = refineEdgeWithSteerableFilter(mask0,image,...
%    doPlot,closureRadius)
%
%INPUT  mask0        : Original mask to be refined.
%       image        : Image to be segmented.
%       doPlot       : 1 to plot masks in the end, 0 otherwise. Refined
%                      masks shown in green, original masks shown in blue.
%
%OUTPUT mask         : Mask (1 inside cell, 0 outside).
%
%Khuloud Jaqaman, November 2011

%% Input

if nargin < 2
    error('refineEdgeWithSteerableFilter: Wrong number of input arguments');
end

if nargin < 3 || isempty(doPlot)
    doPlot = 0;
end

if nargin < 4 || isempty(closureRadius)
    closureRadius = 5;
end

%% Pre-processing

%run steerable filter to enhance edges
[res,theta] = steerableDetector(image, 3, 1.5);

%run Gaussian filter to smoothen noise
imageF = filterGauss2D(image,1.5);

%divide gradient by intensity value to enhance the real cell edge
resOverImage = res./imageF;

%get the non-maximum suppression image to get edge candidates
nmsResOverImage = nonMaximumSuppression(resOverImage,theta);

%extract some connected componenet properties in the nms image
nmsBW = (nmsResOverImage > 0);
[nmsL,numL] = bwlabel(nmsBW,4);
[candLength,candMeanGradient,candVarTheta] = deal(NaN(numL,1));
[nmsMeanGrad,nmsVarTheta] = deal(zeros(size(nmsResOverImage)));
for iL = 1 : numL
    %pixels
    pixL = find(nmsL==iL);
    %properties
    candLength(iL) = length(pixL);
    candMeanGradient(iL) = mean(nmsResOverImage(pixL));
    candVarTheta(iL) = var(theta(pixL));
    %mean gradiant and variance of theta "images"
    nmsMeanGrad(pixL) = candMeanGradient(iL);
    nmsVarTheta(pixL) = candVarTheta(iL);
end

%determine the upper threshold for edge segmentation
cutThreshGrad = prctile(candMeanGradient,95);
% cutThreshTheta = prctile(candVarTheta,85);
cutThreshHigh = cutThreshGrad;

%define range of lower threshold for edge segmentation
cutThreshLow = (1:-0.1:0.3)'*cutThreshGrad;

%get the edges at image boundary of original mask
rangeVal = 25;
maskTmp = bwdist(mask0) <= rangeVal;
mask2 = zeros(size(mask0)+2);
mask2(2:end-1,2:end-1) = maskTmp;
mask2Bound = mask2-imerode(mask2,strel('square',3));
mask0Bound = mask2Bound(2:end-1,2:end-1);
mask0Bound(2:end-1,2:end-1) = 0;

%% Thresholding + post-processing

segmentOK = 0;
jIter = 1;
while ~segmentOK && jIter <= length(cutThreshLow)
    
    %apply gradient threshold to get edge mask
    edgeMask = hysteresisThreshold(nmsMeanGrad,cutThreshHigh,cutThreshLow(jIter));
    
    %     %also apply angle threshold
    %     edgeMask2 = (nmsVarTheta > 0) & (nmsVarTheta < cutThreshTheta);
    %     edgeMask = edgeMask & edgeMask2;
        
    %bridge gaps in edge mask
    edgeMask = bwmorph(edgeMask,'bridge');
    
    %remove isolated edges whose area is not larger than 3 pixels
    regionStats = regionprops(edgeMask,'Area');
    regionArea = vertcat(regionStats.Area);
    indxBad = find(regionArea <= 3);
    L = bwlabel(edgeMask);
    for iBad = indxBad'
        edgeMask(L==iBad) = 0;
    end
    
    %close gaps by image closing
    SE = strel('disk',closureRadius,0);
    edgeMask = imclose(edgeMask,SE);
    
    %skeletonize
    edgeMask = bwmorph(edgeMask,'skel',Inf);
    
    %make the mask by filling holes
    %first add image boundary-related edge
    edgeMask = edgeMask | mask0Bound;
    mask = imfill(edgeMask,'holes');
    
    %get rid of spikes in mask by doing an opening
    SE = strel('disk',2,0);
    mask = imopen(mask,SE);
    
    %if there is > 1 connected component, keep only the largest one
    regionStats = regionprops(mask,'Area');
    regionArea = vertcat(regionStats.Area);
    maxArea = max(regionArea);
    indxBad = find(regionArea<maxArea);
    L = bwlabel(mask);
    for iBad = indxBad'
        mask(L==iBad) = 0;
    end
    
    %check whether segmentation is OK
    if maxArea > 10000
        segmentOK = 1;
    else
        jIter = jIter + 1;
    end

end

%output
if segmentOK
    cutThreshHL = [cutThreshHigh cutThreshLow(jIter)];
else
    mask = zeros(size(mask0));
    cutThreshHL = [cutThreshHigh -1];
end

%% Display

%show original mask and modified mask
if doPlot
    SE = strel('square',3);
    maskEdge0 = mask0 - imerode(mask0,SE);
    maskEdge = mask - imerode(mask,SE);
    imageNorm = (image-min(image(:))) / (max(image(:))-min(image(:)));
    image3 = imageNorm;
    image3(:,:,2) = maskEdge;
    image3(:,:,3) = maskEdge0;
    imtool(image3,[]);
end

%% ~~~ the end ~~~
