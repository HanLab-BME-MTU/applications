function [imgOut, mask, vI, fvI, bfvI] = segCellLCH(img, varargin)
% simple script to pre-process/segment LCH cells for Deep learning.
% use example: [mask vI fvI bfvI]=segCellLCH(imread('./14-May-2017_atcc_s06_t120_x998_y1586_t130_f6.png'));
%
%  INPUT: image
%           'align' : orient the mask along major axis
%           'preview': plot results for debugging
%  OUTPUT: binary mask
%    (optional outputs include image processing steps for debugging)
%
%   
%
% by Andrew R. Jamieson, Oct 2017

ip = inputParser;
ip.addRequired('img', @isnumeric);
ip.addOptional('align',false, @islogical);
ip.addOptional('preview', true, @islogical);
ip.parse(img,varargin{:});
p = ip.Results;

% Core image processing steps
%%%%%%%%%%%%%%%%%%%%%%%
I = mat2gray(img);
gI = imfilter(I, fspecial('gaussian', 5,1));
vI = stdfilt(gI);
fvI = imfilter(vI, fspecial('gaussian', 7,2));
bfvI = imbinarize(fvI, .02);
bfvI = imdilate(bfvI, strel('disk', 2));
maskAll = imclose(bfvI, strel('disk', 7));
maskAll = bwfill(maskAll,'holes');
%%%%%%%%%%%%%%%%%%%%%%%

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

imgOut=mask.*I;

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