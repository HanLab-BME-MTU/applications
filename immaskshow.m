function immaskshow( im , inmasks , displayrange, maskColorMap , maskAlpha )
% immaskshow( im , masks , intrange, maskColorMap , maskAlpha )
%
% Shows masks overlaid on an image.
% NOTE: A pixel can be TRUE in more than one mask.

if isempty(inmasks)
    masks = [];
elseif(isstruct(inmasks))
    masks = inmasks;
elseif iscell(inmasks)
    for i = 1:numel( inmasks )
        masks(i).im = inmasks{i};
    end	   
elseif length(size(im))==length(size(inmasks))
    if length(unique(size(im)==size(inmasks)))==1
        masks(1).im = inmasks;
    end
else    
    error('ERROR: check usage of masks');
end
clear inmasks;
numMasks  = size(masks,2);

if nargin < 5
    maskAlpha(1:numMasks) = 0.5;
end
if nargin < 4
    if(numMasks < 7)
        cMap(1,:) = [1 0 0];
        cMap(2,:) = [0 1 0];
        cMap(3,:) = [0 0 1];
        cMap(4,:) = [1 0 1];
        cMap(5,:) = [0 1 1];
        cMap(6,:) = [1 1 0];
        maskColorMap = cMap(1:numMasks,:);
    end
end
if nargin < 3 || numel(displayrange) ~= 2
    displayrange = double([ min(im(:)) max(im(:)) ]);
end

numColors = size( maskColorMap , 1 );
numAlpha  = length( maskAlpha );
if length( unique( [ numMasks numColors numAlpha ] ) ) ~= 1
   error('error: number of mask images, alpha values and colors are not the same');
end

% Transform image to rgb
maxVal = 255;
volSize = size(im);
imr = zeros(volSize, 'uint8');
img = zeros(volSize, 'uint8');
imb = zeros(volSize, 'uint8');

% To rgb
imr(:,:) = im2uint8(mat2gray(im, displayrange));
img(:,:) = imr(:,:);
imb(:,:) = imr(:,:);

% Add Masks
% fprintf(1, 'Compute RGB Masks ...\n');
for i = 1:numMasks
    masks(i).im = logical(masks(i).im);

    % RGB
    imr(masks(i).im) = (1 - maskAlpha(i)) * imr(masks(i).im) + maxVal * maskAlpha(i) * maskColorMap(i,1);
    img(masks(i).im) = (1 - maskAlpha(i)) * img(masks(i).im) + maxVal * maskAlpha(i) * maskColorMap(i,2);
    imb(masks(i).im) = (1 - maskAlpha(i)) * imb(masks(i).im) + maxVal * maskAlpha(i) * maskColorMap(i,3);
end

% RGB Image
imrgb = cat(3, imr, img, imb);
clear imr img imb masks im;
imshow(imrgb);

% imshow( im , displayrange );
% hold on;
% 
% imsize = size( masks(1).im );
% 
% for i = 1:numMasks
%     imRGBCurMask = zeros( [ imsize 3 ] );
%     for c = 1:3
%         imRGBCurMask( : , : , c ) = masks(i).im * maskColorMap( i , c );
%     end
%     imCurAlpha = masks(i).im * maskAlpha( i );
%     image( imRGBCurMask , 'AlphaData' , imCurAlpha );
%     hold on;
% end
% hold off;
% drawnow;
end