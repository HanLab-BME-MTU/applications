function [] = segmentAndSaveMovie(type, nsize, numf);
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


%% load image files that are to be segmented and make directory for results

ordir = cd;

[oriImageName, oriImagePath] = uigetfile('.tif','Please select first original image for segmentation'); 
oriImageName = strcat(oriImagePath, oriImageName);
oriImageStackList = getFileStackNames(oriImageName);
total_frame_num = length(oriImageStackList);
if nargin>2
    if numf<total_frame_num
        total_frame_num = numf;
        oriImageStackList = oriImageStackList(1:numf);
    end
end

% check whether the \SubregionsMask\ subdirectory exists or not 
[success, msg, msgID] = mkdir(oriImagePath,'SubregionsMask');

if (success ~= 1)
    error(msgID, msg); 
elseif (~isempty(msg))
    fprintf('Subdirectory "SubregionsMask" already exists.\n');
else
    fprintf('Subdirectory "SubregionsMask" has been created.\n');
end

cd(oriImagePath);
cd('SubregionsMask');


% load and display first image
testName = char(oriImageStackList(1));
testImage = imread(testName);

%figure;
imshow(testImage,[]);

[isx,isy] = size(testImage);
objectsize = 0.5*min(isx,isy);
if nargin>1
    if ~isempty(nsize)
        objectsize = nsize;
    end
end

%% if type = cell, look for and load detection data
if strcmp(type,'cell');
    % look for detection data
    cd(oriImagePath);
    detection = [];
    if exist('DetectionStructures')==7
        cd('DetectionStructures')
        open detection.mat;
        det = ans.detection;
    else
        error('no detection.mat file exists for this condition');
    end
    cd(oriImagePath);
    cd('SubregionsMask');
    
    % calculate 'basis mask', which is the overlay of all detected objects
    % first, set central pixels of detected objects to 1
    mask_detection = 0*testImage;
    for t=1:length(det)
        detxpos = det(t).xCoord;
        detypos = det(t).yCoord;
        for n=1:length(detypos)
            mask_detection(round(detypos(n)),round(detxpos(n))) = 1;
        end
    end
    % second, dilate the single pixels to little discs
    se5 = strel('disk',5);
    detmask_detectdilate = imdilate(logical(mask_detection),se5);
    % third, fill the resulting area and filter with Gaussian profile to
    % create a continuous area
    detmask_filled = imfill(detmask_detectdilate,'holes');
    detmask_Gfiltered = Gauss2DBorder(detmask_filled,10);
end


%% loop over all images and perform segmentation

% loop over all images
for i=1:total_frame_num
    
    % load current image
    tempname = char(oriImageStackList(i));
    currImage = imread(tempname);
    [cisx,cisy,cisz] = size(currImage);
    if cisz>1
        currImage = currImage(:,:,1);
    end
    
    % if segmentation is to be optimized for patterned surface
    if strcmp(type,'pattern')
        mask = segmentArea_pattern(currImage, objectsize);
    elseif strcmp(type,'cell');
        % look for detection data
        mask = segmentArea_cell(currImage, objectsize, detmask_Gfiltered);
    end

    currentname = sprintf('mask%04d',i);
    cfilename = [currentname,'.tif'];
       
    imwrite(mask,cfilename,'tif');
      
end

cd(ordir);

end % of function