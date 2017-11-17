function [imgOut, mask, vI, fvI, bfvI] = segCellLCH(img, varargin)
% simple script to pre-process/segment LCH cells for Deep learning.
% use example: [mask vI fvI bfvI]=segCellLCH(imread('./14-May-2017_atcc_s06_t120_x998_y1586_t130_f6.png'));
%
%  INPUT: image
%           'align' : orient the mask along major axis
%           'preview': plot results for debugging
%             'avgBG':    default is false for using zeros, 
%                         otherwise fills the non-mask with the avg of the 
%                         original image background
%             'varFilterOut': outputs the variance filter of the segmentated 
%                         image region, note variance filter is done first on entire 
%                         image the mask is cut from this
%             'centerMass": place object in the center of the image and pad
%             with zero/background
%             
%  OUTPUT: binary mask
%    (optional outputs include image processing steps for debugging)
%
%   
%
% by Andrew R. Jamieson, Oct 2017

ip = inputParser;
ip.addRequired('img', @isnumeric);
ip.addOptional('align',false, @islogical);
ip.addOptional('preview', false, @islogical);
ip.addOptional('avgBG', false, @islogical);
ip.addOptional('varFilterOut', false, @islogical);
ip.addOptional('centerMass', false, @islogical);
ip.parse(img,varargin{:});
p = ip.Results;

% Core image processing steps
%%%%%%%%%%%%%%%%%%%%%%%
I = mat2gray(img);
gI = imfilter(I, fspecial('gaussian', 5,.25));
vI = stdfilt(gI);
fvI = imfilter(vI, fspecial('gaussian', 5,1));
bfvI = imbinarize(fvI, .02);
bfvI = imdilate(bfvI, strel('disk', 5));
maskAll = imclose(bfvI, strel('disk', 5));
% maskAll = imerode(maskAll, strel('disk', 4));
% maskAll = imerode(maskAll, strel('disk', 2));
maskAll = bwfill(maskAll,'holes');
%%%%%%%%%%%%%%%%%%%%%%%

% I = mat2gray(img);
% gI = imfilter(I, fspecial('gaussian', 5,1));
% vI = stdfilt(gI);
% fvI = imfilter(vI, fspecial('gaussian', 7,2));
% bfvI = imbinarize(fvI, .02);
% bfvI = imdilate(bfvI, strel('disk', 2));
% maskAll = imclose(bfvI, strel('disk', 7));
% maskAll = bwfill(maskAll,'holes');

% See how many objects
CC = bwconncomp(maskAll);

% check if in center
% if yes, keep just this object
% if multiple, select one that overlaps with center
if CC.NumObjects > 1
    centerPts = round(size(maskAll)./2);
    mask = bwselect(maskAll,centerPts,centerPts);

    % check if mask present at center.
    if isempty(find(mask,1))
        % if not, just take the larget
        mask = bwareafilt(maskAll,1);
    end
else
    mask = maskAll;
end


% check if 97% covering image, then
% rp = regionprops(mask);
% if rp.Area/size(maskAll,1)^2  >= .9
%     mask = 0;
% end
if p.varFilterOut
%     Ivar = rangefilt(I);
    Ivar = stdfilt(I);
    imgFG = mask.*Ivar;    
    
else
    imgFG = mask.*I;    
end


if p.avgBG && ~p.varFilterOut % (note, just use zeros for var filter outputs
    rp = regionprops(mask);
    if rp.Area ~= size(maskAll,1)^2
        imgBG = ~mask.*I;
        imgBG = ~mask .* mean(imgBG(~mask));
        imgOut = imgFG + imgBG;
    else
        imgOut = imgFG;
    end
    
else
    imgOut = imgFG;
end

if p.align
    disp('re-orienting mask')
    rp = regionprops(mask,'orientation');
    rp.Orientation
    imgOut = imrotate(imgOut,-1*rp.Orientation);
end

if p.centerMass
    rp = regionprops(mask);
    
    [r c] = size(mask);
    rShift = round(r/2 - rp.Centroid(2));
    cShift = round(c/2 - rp.Centroid(1));

    % Call circshift to move region to the center.
    imgOut = circshift(imgOut, [rShift cShift]);
    mask = circshift(mask, [rShift cShift]);
end


if p.align
    disp('re-orienting mask')
    rp = regionprops(mask,'orientation');
    rp.Orientation
    imgOut = imrotate(imgOut,-1*rp.Orientation);
end


if p.preview 
    % preview results
    figure; imshow(imgOut);
    figure;
    subplot(3,2,1);imshow(bfvI,[]); title('after thre & dilation');
    subplot(3,2,2);imshow(vI,[]);title('variance filter');
    subplot(3,2,3);imshow(fvI,[]); title('smoothing');
    subplot(3,2,4);imshow(imgOut,[]);title('final');
    subplot(3,2,5);imshow(maskAll,[]);title('all object mask');
    subplot(3,2,6);imshow(I,[]); title('original');
end

% Works roughly for 128x128 downsampling
% vI = stdfilt(I);
% fvI = imfilter(vI, fspecial('gaussian', 7,3));
% bfvI = imbinarize(fvI, .02);
% bfvI = imdilate(bfvI, strel('disk', 3));
% mask = imclose(bfvI, strel('disk', 7));