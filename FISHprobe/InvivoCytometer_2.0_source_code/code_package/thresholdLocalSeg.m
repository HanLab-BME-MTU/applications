function  [level,bw_out,varargout] = thresholdLocalSeg(imageIn,choice_of_threshold, level_local_radius, pace, lowerbound, varargin)
% local thresholding based on thresholding levels determined by the
% thrsholding method of choice
%
% [level,bw_out] = thresholdLocal(imageIn,choice_of_threshold, level_local_radius, pace, showPlots)
% 
% This function selects a global threshold for the input fluorescence image using
% the thresholding method the user askes for, indicated in choice_of_threshold, 
% then find the local thresholds also
% decided the chosen thresholding method to get a local threshold map. 
% After some smoothing this threshold map is used to segment the input% image.
% 
% Input:
% 
%   imageIn:            2D input image to be thresholded.
%   choice_of_threshold: string input "Otsu","Rosin",or "FluorescenceImage"
%   level_local_radius: the radius of local patch
%   pace:               the pace to calculate local threshold, mostly to
%                       speed up the process which is computational expensive
%   lowerbound:         the percentage of global threshold the user want to
%                       set as the lowerbound of the local thresholding
%   showPlots:          If true, a plot of the histogram and an overlay of the mask
%                       on the image will be shown. 
%
% Output:
% 
%   level - The intensity value selected for global thresholding.
%   bw_out - the output segmentation result
%
% Liya Ding, 12/2012

ip=inputParser;
ip.addRequired('imageIn',@isnumeric);
ip.addRequired('choice_of_threshold', ...
               @(x) ((ischar(x) && ismember(x, {'Otsu', 'Rosin', 'FluorescenceImage'})) || ...
                      isa(x, 'function_handle')));
ip.addRequired('level_local_radius',@isnumeric);
ip.addRequired('pace',@isnumeric);
ip.addRequired('lowerbound',@isnumeric);
ip.addOptional('showPlots',0,@isnumeric)
ip.parse(imageIn, choice_of_threshold, level_local_radius, pace, lowerbound, varargin{:});
showPlots=ip.Results.showPlots;


% level_local_radius = 15;
% pace = 3;
half_pace = round(pace-1)/2;

% Convert to double if necessary
imageIn = double(imageIn);

% Calculate the local and global threshold using the thresholding method
% indicated in choice_of_threshold
[level_img, level_whole] = thresholdLocalCalculate(imageIn,choice_of_threshold,level_local_radius, pace);

% Smooth the threshold map
H = fspecial('gaussian',4*pace+1,half_pace);
level_img = imfilter(level_img,H,'replicate','same')/2;

% Segmentation using the threshold map
imageMask = imageIn >= level_img;

% Get a global mask to eliminate segmentation of detailed noise in the
% background, with a slight lowered threshold to include some boundary
% parts, fill holes and dilate a little to avoid masking off target region

% second_group =imageIn(find(imageIn>level_whole));
% new_level_for_second_group = level_whole+lowerbound/100*std(second_group);
% imageMask_whole = imageIn>new_level_for_second_group;

imageMask_whole = imageIn >= level_whole*lowerbound/100;
% imageMask_whole = imfill(imageMask_whole,'holes');
% imageMask_whole = imdilate(imageMask_whole,ones(5,5));

% The final segmentation is the intersect of both mask.
bw_out =  imageMask.*imageMask_whole;

% The output level is the gobal threshold 
level = level_whole;

if(showPlots==1)
    figure(1); hold off;
    imagesc(imageIn); colormap(gray); hold on
    contour(bw_out,'r')
end

if nargout > 2
    minThresh = level_whole*lowerbound/100;
    level_img( level_img < minThresh ) = minThresh;
    varargout{1} = level_img;
end

