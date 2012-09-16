function [mask,nmsThresh] = refineEdgeWithSteerableFilter(mask0,image,...
    doPlot,nmsThresh,highLowFactor,closureRadius)
%REFINEEDGEWITHSTEERABLEFILTER refines cell edge using intensity gradients obtained from a steerable line filter
%
%SYNOPSIS [mask,nmsThresh] = refineEdgeWithSteerableFilter(mask0,image,...
%    doPlot,nmsThresh,highLowFactor,closureRadius)
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

if nargin < 4 || isempty(nmsThresh)
    nmsThresh = [];
end

if nargin < 5 || isempty(highLowFactor)
    highLowFactor = [1.2 0.6];
end

if nargin < 6 || isempty(closureRadius)
    closureRadius = 15;
end

%% Range of image pixels to consider

rangeVal = 15;
D1 = ~mask0 & (bwdist(mask0) <= rangeVal);
D2 = mask0 & (bwdist(~mask0) <= rangeVal);
maskInRange = D1 | D2;

%% Filtering and thresholding

%run steerable filter to enhance edges
[res,theta] = steerableDetector(image, 3, 1.5);

%divide gradient by intensity value to enhance the real cell edge
resOverImage = res./image;

%do non-maximum suppression on the gradient-to-intensity ratio
nms = nonMaximumSuppression(resOverImage,theta);

%keep only value in the correct image range
nms(~maskInRange) = 0;

%determine thresholds using the Rosin algorithm
nmsVec = nms(nms~=0);
if isempty(nmsThresh)
    [~,nmsThresh] = cutFirstHistMode(nmsVec,0);
    %     nmsThresh = max(nmsThresh,prctile(nmsVec,90));
end
nmsThreshHigh = highLowFactor(1) * nmsThresh;
nmsThreshLow  = highLowFactor(2) * nmsThresh;

%add the edges at image boundary of original mask to the nms image
maskTmp = bwdist(mask0) <= rangeVal;
mask2 = zeros(size(mask0)+2);
mask2(2:end-1,2:end-1) = maskTmp;
mask2Bound = mask2-imerode(mask2,strel('square',3));
mask0Bound = mask2Bound(2:end-1,2:end-1);
mask0Bound(2:end-1,2:end-1) = 0;
mask0Bound = mask0Bound*nmsThreshHigh*1.1;
nms = nms + mask0Bound;

%apply threshold to nms image to get cell edge
edgeMask = hysteresisThreshold(nms,nmsThreshHigh,nmsThreshLow);

%% Morphological operations to get a continuous edge without spikes

%bridge gaps in edge mask and skeletonize
edgeMask = bwmorph(edgeMask,'bridge');
edgeMask = bwmorph(edgeMask,'skel',Inf);

%remove isolated edges whose area is not larger than 3 pixels
regionStats = regionprops(edgeMask,'Area');
regionArea = vertcat(regionStats.Area);
indxBad = find(regionArea <= 3);
L = bwlabel(edgeMask);
for iBad = indxBad'
    edgeMask(L==iBad) = 0;
end

%do a closure operation in case there are gaps
SE = strel('disk',closureRadius,0);
% edgeMaskEnds = bwmorph(edgeMask,'endpoints') | mask0Bound;
% edgeMaskEnds = imclose(edgeMaskEnds,SE);
% edgeMask = edgeMask | edgeMaskEnds;
edgeMask = imclose(edgeMask,SE);

%fill holes then remove insides in order to thin the edge
edgeMaskTmp = zeros(size(edgeMask));
edgeMaskTmp(2:end-1,2:end-1) = edgeMask(2:end-1,2:end-1);
edgeMaskTmp = imfill(edgeMaskTmp,'holes');
edgeMask = edgeMaskTmp | mask0Bound;
edgeMask= bwmorph(edgeMask,'remove');

%fill holes manually
imageNorm = (image-min(image(:))) / (max(image(:))-min(image(:)));
SE = strel('square',3);
maskEdge0 = mask0 - imerode(mask0,SE);
figure
imshow(cat(3,imageNorm,edgeMask,maskEdge0));
figure
mask = imfill(edgeMask);
close all

%get rid of spikes in mask by doing an opening
SE = strel('disk',2);
mask = imopen(mask,SE);

%add the image boundary-related edge again and fill holes again to take
%care of corners
%get rid of spikes again
mask = mask | mask0Bound;
mask = imfill(mask,'holes');
SE = strel('disk',2);
mask = imopen(mask,SE);

%if there is > 1 connected component, keep only the largest one
regionStats = regionprops(mask,'Area');
regionArea = vertcat(regionStats.Area);
indxBad = find(regionArea<max(regionArea));
L = bwlabel(mask);
for iBad = indxBad'
    mask(L==iBad) = 0;
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
