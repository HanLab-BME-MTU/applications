function [mask] = segmentArea_pattern(image,nsize);
% segment patch(es) in inhomogeneously illuminated image using adaptive
% thresholding
% SYNOPSIS   [mask] = segmentArea(image,nsize);
%       
% INPUT     image:    grayscale image
%           nsize:    window size for adaptive thresholding; this value
%                     should be chosen such that a window of this size
%                     always contains part of an object/patch (as opposed
%                     to only containing background area), so the optimum
%                     window size depends on patrch size and density
% OUTPUT    mask:     binary mask image
% 
%    

[sx, sy] = size(image);

% initialize results
mask_total = segmentSigBG(image);
mask_adaptive = 0*mask_total;

% strategy:
% calculate threshold in different segments of the image
% then interpolate the threshold values across the entire image
% then threshold the mask

% initialize grid points using the specified window size; grid points are
% the center of the local sampling window which is used to calculate the
% local threshold value
[px,py] = meshgrid(1:nsize:sx+nsize, 1:nsize:sy+nsize); 
% grid of local threshold values
pz = nan*px;

% loop over all grid points
for i=1:length(px(:))
    
    lx = px(i);
    ly = py(i);
    % window coordinates for intensity scaling - take care for cutoff by
    % the edges of the image
    xmins = max( (lx-nsize), 1);
    xmaxs = min( (lx+nsize), sx);
    ymins = max( (ly-nsize), 1);
    ymaxs = min( (ly+nsize), sy);
    
    % partial image cropped by this window
    image_partial = image(xmins:xmaxs,ymins:ymaxs);
    % partial image pixel coordinates
    [xip,yip] = ndgrid(xmins:xmaxs,ymins:ymaxs);
    % list of x,y-coordinates
    Wmat = [xip(:) yip(:)];
    % x,y-coordinates of central grid point
    Pmat = [lx ly]';
    % distances of all crop window points from central point
    distvec = dist(Wmat,Pmat);
    % restrict analysis to points within nsize distance
    usepos = find(distvec <= nsize);
    % determine threshold for the specified subset of points
    [mask_partial,threshold] = segmentSigBG(image_partial,usepos);
    % enter the local threshold into threshold grid    
    pz(i) = threshold;

end

% surf(px,py,pz);
% construct new meshgrid matching the dimensions of the image exactly
[pxi,pyi] = meshgrid(1:sx, 1:sy); 
% interpolate the adaptive threshold values to all positions in the new
% meshgrid
pzi = interp2(px,py,pz,pxi,pyi);
% threshold input image with this interpolated adaptive value
mask_adaptive = image > pzi';


%figure;


%% ==============================================================
% now convert thresholded features into continuous closed patches

se1 = strel('disk',1);
se5 = strel('disk',5);
se10 = strel('disk',10);
se20 = strel('disk',20);

% fill holes inside the patch with 'fill' or 'close'
mask_filled = imfill(double(mask_adaptive),'holes');

% remove single noisy background pixels outside of patches with 'open'
mask_opened = imopen(double(mask_filled),se1);

% dilate and erode to fill remaining holes and smooth the edge
mask_dil = imdilate(double(mask_opened),se10);
mask_erode = imerode(double(mask_dil),se10);

mask = mask_erode;


% mask overlay (ONLY for display purposes)
imin = min(image(:));
imax = max(image(:));
mask_fin = 0.5*(1-mask)+( double(image-imin)/double(imax-imin) );


subplot(2,4,1); imshow(image,[]);       title('original');
subplot(2,4,2); imshow(pzi',[]);        title('adapt. threshold');
subplot(2,4,3); imshow(mask_total,[]);  title('total thresholded');
subplot(2,4,4); imshow(mask_adaptive,[]); title('adapt.thresholded');
subplot(2,4,5); imshow(mask_filled,[]); title('filled');
subplot(2,4,6); imshow(mask_opened,[]); title('opened');
subplot(2,4,7); imshow(mask_erode,[]);    title('dilated+eroded');
subplot(2,4,8); imshow(mask_fin,[]);    title('overlay');


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