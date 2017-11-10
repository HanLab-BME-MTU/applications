function [mask] = segmentArea_cell(image,nsize,mask_det);
% segment cell area in the image
% SYNOPSIS   [mask] = segmentArea(image,nsize);
%       
% INPUT     image:    grayscale image
%           nsize:    approximate size of patch features (in pixels)
%           mask_det: mask from overlaid detected object positions  
% OUTPUT    mask:     binary mask image
% 
%    

[sx, sy] = size(image);

% threshold original image
[mask_total] = segmentSigBG(image);

% fill, filter and scale detection image
image_detscaled = double(max(image(:))) * mask_det;
image_add = double(image) + image_detscaled;

% threshold combined image
[mask_detadd] = segmentSigBG(image_add);

% combine masks
mask_combine = mask_detadd | mask_total;


%% ==============================================================
% now convert thresholded features into continuous closed patches

% fill holes inside the patch with 'fill' or 'close'
se5 = strel('disk',5);
mask_filled = imfill(logical(mask_combine),'holes');

% remove single noisy background pixels outside of patches with 'open'
se1 = strel('disk',1);
mask_opened = imopen(logical(mask_filled),se1);

% dilate to fill remaining holes and smooth the edge
se5 = strel('disk',5);
mask_dil = imdilate(logical(mask_opened),se5);

mask = mask_dil;

subplot(2,3,1); imshow(image,[]);       title('original');
subplot(2,3,2); imshow(mask_detadd,[]); title('det add');
subplot(2,3,3); imshow(mask_combine,[]);title('combine');
subplot(2,3,4); imshow(mask_filled,[]);title('filled');
subplot(2,3,5); imshow(mask_opened,[]);    title('opened');
subplot(2,3,6); imshow(mask_dil,[]); title('dilated');

pause(0.1);


end

    
    
%% subfunctions

function [mask, threshold] = segmentSigBG(image,usepos);

% segment intensities into two levels - features vs background using kmeans
% clustering
intvec = double(image(:));
if nargin>1
    intvec = double(image(usepos));
end

[DX,cloc] = kmeans(intvec,2,'emptyaction','singleton');

% positions of levels 1 and 2
bool1 = DX-1;
bool2 = ~bool1;

f1 = find(bool1);
f2 = find(bool2);

% threshold that separates the two levels (not knowing a priori which of
% the two levels is feature, and which is background)
imin1 = min(intvec(f1));
imin2 = min(intvec(f2));

threshold = max(imin1,imin2);
absmin = min(imin1,imin2);

% thresholded image
mask = 0*image;
mask = (image>=threshold);

end % of subfunction