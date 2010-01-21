function [] = segmentAndSaveMovie(type, nsize, numf)
% the function uploads tiff files, performs segmentation, and the saves the
% segmented mask files into a folder called Subregionsmask
% SYNOPSIS   [] = segmentAndSaveMovie();
%       
% INPUT     type = string, can be
%           'pattern'  = segmentation is optimized for patterned surface
%                    (regular patches of similar size but with possible
%                     uneven intensity)
%           'cell'     = segmentation is optimized for cell arae (very uneven
%                     intensities, possible clathrin pits)
%           nsize = (OPT) sample size for adaptive threshold; if none is
%                   provided, it's automatically set to half of the smaller
%                   window side length
%           numf =  (OPT) number of frames; if none is provided, it's
%                   automatically  set to all available
%
% Last Modified: Francois Aguet, 01/21/2010


% load image files that are to be segmented and make directory for results

[imageName, imagePath] = uigetfile('.tif', 'Please select first original image for segmentation'); 
imageName = strcat(imagePath, imageName);
oriImageStackList = getFileStackNames(imageName);
total_frame_num = length(oriImageStackList);
if nargin>2 && numf<total_frame_num
    total_frame_num = numf;
    oriImageStackList = oriImageStackList(1:numf);
end

% check whether the \SubregionsMask\ subdirectory exists or not 
[success, msg, msgID] = mkdir(imagePath,'SubregionsMask');

if (success ~= 1)
    error(msgID, msg); 
elseif (~isempty(msg))
    fprintf('Subdirectory "SubregionsMask" already exists.\n');
else
    fprintf('Subdirectory "SubregionsMask" has been created.\n');
end

% load and display first image
testImage = imread([imagePath 'SubregionsMask' filesep char(oriImageStackList(1))]);

%figure;
imshow(testImage,[]);

[isx,isy] = size(testImage);
objectsize = 0.5*min(isx,isy);
if nargin>1 && ~isempty(nsize)
    objectsize = nsize;
end

% if type = cell, look for and load detection data
if strcmp(type, 'cell');
    % look for detection data
    detectionPath = [imagePath 'DetectionStructures' filesep 'detection.mat'];
    if (exist('detectionPath', 'file')==2)
        load(detectionPath);
    else
        error('no detection.mat file exists for this condition');
    end
    
    % calculate 'basis mask', which is the overlay of all detected objects
    % first, set central pixels of detected objects to 1
    mask_detection = zeros(size(testImage));
    for t = 1:length(detection)
        for n = 1:length(detection(t).xCoord)
            mask_detection(round(detection(t).yCoord(n)),round(detection(t).xCoord(n))) = 1;
        end
    end
    % second, dilate the single pixels to little discs
    detmask_detectdilate = imdilate(logical(mask_detection), strel('disk',5));
    % third, fill the resulting area and filter with Gaussian profile to
    % create a continuous area
    detmask_filled = imfill(detmask_detectdilate, 'holes');
    detmask_Gfiltered = Gauss2D(detmask_filled, 10);
end

% loop over all images and perform segmentation
for i=1:total_frame_num
    
    % load current image
    currImage = imread([imagePath 'SubregionsMask' filesep char(oriImageStackList(i))]);
    if size(currImage, 3)>1
        currImage = currImage(:,:,1);
    end
    
    % if segmentation is to be optimized for patterned surface
    if strcmp(type,'pattern')
        mask = segmentArea_pattern(currImage, objectsize);
    elseif strcmp(type,'cell');
        mask = segmentArea_cell(currImage, objectsize, detmask_Gfiltered);
    end
    currentname = sprintf('mask%04d',i);
    imwrite(mask, [currentname,'.tif'], 'tif');      
end