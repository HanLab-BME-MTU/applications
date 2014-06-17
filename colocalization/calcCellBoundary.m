function [maskList,mask] = calcCellBoundary(image,doPlot)
%CALCCELLBOUNDARY applies mutli-otsu threshold to determine cell boundary
% Mask is best suited for images in which cell is large relative to image
% dimensions
% Synopsis: [maskList,mask] = calcCellBoundary(image)
% Input:
%   image - 
%   doPlot: 1 to show result, 0 otherwise. Default: 0.
%
% Output:
%   maskList- n x 2 list of all pixels within cellBoundary
%
%   mask - binary

if nargin < 2 || isempty(doPlot)
    doPlot = 0;
end

%keep original image for plotting
image0 = image;

% Apply wide gaussian filter
 gaussF = fspecial('gaussian',[20 20],10); 
 image = imfilter(image,gaussF,'same');

%Apply multi-Otsu threshold on image
thresh = multithresh(image,4);
test = image >= thresh(1);
test = imfill(test,'holes');


% Find Connected Components and their areas
[L, ~] = bwlabel(test, 8);
STATS = regionprops(L, 'Area');


%Find all areas and then only get indices of components that have at least
%some fraction of the largest area (default: 0.5)
areaInfo = cat(1,STATS.Area);
sThreshold = find(areaInfo >= (max(areaInfo)/2));


% Create new matrix with only the selected components, there's probably an
% easier way of doing this...
for i =1:length(sThreshold)
    mask(:,:,i) = (L ==sThreshold(i));
    mask(:,:,1) = mask(:,:,1)+mask(:,:,i);
end

[i,j] = find(mask(:,:,1));
maskList(:,1) = i;
maskList(:,2) = j;

if doPlot
    figure
    imshow(image0,[]);
    hold on;
    maskBounds = bwboundaries(mask);
    cellfun(@(x)(plot(x(:,2),x(:,1),'r','LineWidth',1)),maskBounds);
end

end