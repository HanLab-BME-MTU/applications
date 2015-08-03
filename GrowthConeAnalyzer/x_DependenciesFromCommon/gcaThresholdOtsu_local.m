function [level,bw_out] = gcaThresholdOtsu_local(imageIn,level_local_radius, pace, varargin)
% local thresholding based on global thresholding level using Otsu's method
%
% [level,bw_out] = thresholdOtsu_local(imageIn,level_local_radius, pace, showPlots)
% 
% This function selects a threshold for the input fluorescence image using
% Ostu method(function thresholdOtsu()), then find the local thresholds also
% decided by Ostu method(function thresholdOtsu()) to get a local threshold
% map. After some smoothing this threshold map is used to segment the input
% image.
% 
% Input:
% 
%   imageIn:            2D input image to be thresholded.
%   level_local_radius: the radius of local patch
%   pace:               the pace to calculate local threshold, mostly to
%                       speed up the process which is computational expensive
%   showPlots:          If true, a plot of the histogram and an overlay of the mask
%                       on the image will be shown. 
%
% Output:
% 
%   level - The intensity value selected for global thresholding.
%
% Liya Ding, 1/2012

ip=inputParser;
ip.addRequired('imageIn',@isnumeric);
ip.addRequired('level_local_radius',@isnumeric);
ip.addRequired('pace',@isnumeric);
ip.addOptional('showPlots',0,@isnumeric)
ip.parse(imageIn, level_local_radius, pace, varargin{:});
showPlots=ip.Results.showPlots;

% level_local_radius = 15;
% pace = 3;
half_pace = round(pace-1)/2;

% Convert to double if necessary
imageIn = double(imageIn);

% A light gaussian smoothing
H1 = fspecial('gaussian',3);
imageIn =  imfilter(imageIn, H1,'same','replicate');

% Find nonzero values (due to masking)
nzInd = find(imageIn);

% Get minumum and maximum pixel values in image
minSignal = min(imageIn(nzInd));
maxSignal = max(imageIn(nzInd));

% Normalize nonzero value between 0 and 1
imageInNorm = zeros(size(imageIn));
imageInNorm(nzInd) = (imageIn(nzInd)- minSignal) / (maxSignal - minSignal);

% Gaussian smoothing of the normalized image, to get a better global mask
H = fspecial('gaussian',7);
sm_imageInNorm = imfilter(imageInNorm, H,'same','replicate');
level = graythresh(sm_imageInNorm);
level_sm_norm = level;

% the glocal level determined for input image
level_whole = level*(maxSignal - minSignal)+minSignal;

% Initialized the threshold map(level_img) as setting everywhere as the global level
level_img = imageIn*0+level_whole;

% For every patch according to the pace(grid)
for img_x = level_local_radius + 1 : pace : size(level_img,1) - level_local_radius
    for img_y = level_local_radius + 1 : pace : size(level_img,2) - level_local_radius

        % Cropping a local image path
        local_img = imageIn(img_x - level_local_radius : img_x + level_local_radius ,img_y - level_local_radius:img_y + level_local_radius );
        
        % Find nonzero values (due to masking)
        nzInd = find(local_img);
        
        % Get minumum and maximum pixel values in image
        minSignal = min(local_img(nzInd));
        maxSignal = max(local_img(nzInd));
        
        % Normalize nonzero value between 0 and 1
        imageInNorm = zeros(size(local_img));
        imageInNorm(nzInd) = (local_img(nzInd)- minSignal) / (maxSignal - minSignal);
        
        % Find the threshold for this patch
        level = graythresh(imageInNorm);
        level = level*(maxSignal - minSignal)+minSignal;
        
        % level could be empty if input patch is all background zeros
        if(isempty(level))
            level_img(img_x-half_pace:img_x+half_pace,img_y-half_pace:img_y+half_pace) = level_whole;
        else
            level_img(img_x-half_pace:img_x+half_pace,img_y-half_pace:img_y+half_pace) = level;
        end
    end
end

% Smooth the threshold map
H = fspecial('gaussian',pace);
level_img = imfilter(level_img,H,'replicate','same');

% Segmentation using the threshold map
imageMask = imageIn >= level_img;

% Get a global mask to eliminate segmentation of detailed noise in the
% background, with a slight lowered threshold to include some boundary
% parts, fill holes and dilate a little to avoid masking off target region
imageMask_whole = sm_imageInNorm >= level_sm_norm*0.85;
imageMask_whole = imfill(imageMask_whole,'holes');
imageMask_whole = imdilate(imageMask_whole,ones(5,5));

% The final segmentation is the intersect of both mask.
bw_out =  imageMask.*imageMask_whole;

% The output level is the gobal threshold 
level = level_whole;

if(showPlots==1)
    figure(1); hold off;
    imagesc(imageIn); colormap(gray); hold on
    contour(bw_out,'r')
end
