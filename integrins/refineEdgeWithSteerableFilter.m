function mask = refineEdgeWithSteerableFilter(mask0,image,doPlot)
%REFINEEDGEWITHSTEERABLEFILTER refines cell edge using intensity gradients obtained from a steerable line filter
%
%SYNOPSIS mask = refineEdgeWithSteerableFilter(mask0,image,doPlot)
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

if nargin ~= 3
    error('refineEdgeWithSteerableFilter: Wrong number of input arguments');
end

%% Filtering and thresholding

%run steerable filter to enhance edges
[res,theta] = steerableDetector(image, 3, 1.5);

%divide gradient by intensity value to enhance the real cell edge
resOverImage = res./image;

%do non-maximum suppression on the gradient-to-intensity ratio
nms = nonMaximumSuppression(resOverImage,theta);

%determine thresholds using the Rosin algorithm
nmsVec = nms(nms~=0);
[~,nmsThreshRosin] = cutFirstHistMode(nmsVec,0);
nmsThreshRosin = max(nmsThreshRosin,prctile(nmsVec,90));
nmsThreshHigh = 1.5 * nmsThreshRosin;
nmsThreshLow  = 0.8 * nmsThreshRosin;

%apply threshold to nms image to get cell edge
edgeMask = hysteresisThreshold(nms,nmsThreshHigh,nmsThreshLow);

%% Morphological operations to get a continuous edge without spikes

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

%determine number of connected components
[~,numEdges] = bwlabel(edgeMask);

%do a closure operation if there is > 1 connected component
if numEdges > 1
    closureRadius = 3;
    SE = strel('disk',closureRadius);
    edgeMask = imclose(edgeMask,SE);
end

%combine the original mask with the one from the edge detection
mask = logical(edgeMask) | logical(mask0);

%fill holes
mask = imfill(mask,'holes');

%get rid of spikes in mask by doing an erosion followed by a dilation
SE = strel('disk',10);
mask = imerode(mask,SE);
mask = imdilate(mask,SE);

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
