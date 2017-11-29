function [maskBlobs,labels,maskBlobsVis] = segmentBlobs_locmax(image,...
    thresholdMethod,methodValue,filterNoise,filterBackground,minSize,plotRes,mask) 
%segmentBlobs_locmax segments blobs in 2Detection_IdxD images via various 
%thresholding methods and use local maxima to catch the dim objects
%
%SYNOPSIS [maskBlobs, imageFilteredMinusBackgroundNorm,labels] = blobSegmentThreshold(image,...
%  threshold,methodValue,filterNoise,filterBackground,minSize,plotRes,mask)
%
%INPUT  image     : 2D image to be segmented.
%       thresholdMethod: 
%                   'otsu' for Otsu.
%                   'rosin' for Rosin.
%                   'minmax' for first minimum after first maximum.
%                   'prct' to use a certain percentile.
%                   'user' to use a threshold input by user.
%                   Optional. Default: 'otsu'.
%       methodValue: Needed only if thresholdMethod = 'prct' or 'user'.
%                    If 'prct', then this is the percentile to use.
%                    Optional. Default: 90.
%                    If 'user', then this is the threshold value to use.
%                    Optional. Default: 0.9.
%       filterNoise: Either 0 to not filter noise or filter sigma > 0 to
%                    filter noise.
%                    Optional. Default: 1.
%       filterBackground: Either 0 to not filter background or filter sigma
%                         > 0 to filter background.
%                         Optional. Default: 10.
%       minSize   : Minimum size of a blob. 
%                   Optional. Default: 20 pixels.
%       plotRes   : 1 to plot segmentation results, 0 otherwise.
%                   Optional. Default: 0.
%       mask      : Binary mask. Optional. If not provided, the whole image
%                   domain is segmented.
%
%OUTPUT maskBlobs : mask including segmentations and local maxima.  
%     
%
%      Labels  : Labels of objects in maskBlobs.
%   
%Khuloud Jaqaman October 2012
%Jesus Vega July 2017
 

%% Input

%check number of input arguments
if nargin < 1
    disp('Please enter at least image to be segmented');
    return
end

%thresholding method
if nargin < 2 || isempty(thresholdMethod)
    thresholdMethod = 'otsu';
end

%method value
if nargin < 3 || isempty(methodValue)
    switch thresholdMethod
        case 'prct'
            methodValue = 90;
        case 'user'
            methodValue = 0.9;
        otherwise
            methodValue = [];
    end
end

%noise filtering
if nargin < 4 || isempty(filterNoise)
    filterNoise = 1;
end

%background filtering
if nargin < 5 || isempty(filterBackground)
    filterBackground = 10;
end

%minimum blob size
if nargin < 6 || isempty(minSize)
    minSize = 20;
end

%plot results
if nargin < 7 || isempty(plotRes)
    plotRes = 0;
end

%mask
if nargin < 8 || isempty(mask)
    mask = ones(size(image));
end

if ~logical(mask)
    error('Mask must be a logical image.');
end
    
%% Segmentation

%make sure that image is in double format
image = double(image);

%remove noise by filtering image with a narrow Gaussian
if filterNoise > 0
    imageFiltered = filterGauss2D(image,filterNoise);
else
    imageFiltered = image;
end 

%estimate background by filtering image with a wide Gaussian
if filterBackground > 0
    imageBackground = filterGauss2D(image,filterBackground);
else
    imageBackground = zeros(size(image));
end

%calculate noise-filtered and background-subtracted image
imageFilteredMinusBackground = imageFiltered - imageBackground;

%crop image
%imageFilteredMinusBackground = imageFilteredMinusBackground .* mask;

%find nonzero values (due to masking)
nzInd = find(imageFilteredMinusBackground);

%get minumum and maximum pixel values in image
minSignal = min(imageFilteredMinusBackground(nzInd));
maxSignal = max(imageFilteredMinusBackground(nzInd));

%normalize nonzero value between 0 and 1
imageFilteredMinusBackgroundNorm = zeros(size(imageFilteredMinusBackground));
imageFilteredMinusBackgroundNorm(nzInd) = (imageFilteredMinusBackground(nzInd) - minSignal) / (maxSignal - minSignal);

%estimate the intensity level to use for thresholding the image
switch thresholdMethod
    case 'otsu'
        try
            level = graythresh(imageFilteredMinusBackgroundNorm(nzInd));
        catch %#ok<CTCH>
            disp(['Method failed to determine a threshold (otsu, noise ' ...
                num2str(filterNoise) ', background ' ...
                num2str(filterBackground) '). Exiting code.'])
            return
        end
    case 'rosin'
        try
            [dummy,level] = cutFirstHistMode(imageFilteredMinusBackgroundNorm(nzInd),0); %#ok<ASGLU>
        catch %#ok<CTCH>
            disp(['Method failed to determine a threshold (rosin, noise ' ...
                num2str(filterNoise) ', background ' ...
                num2str(filterBackground) '). Exiting code.'])
            return
        end
    case 'minmax'
        try
            level = thresholdFluorescenceImage(imageFilteredMinusBackgroundNorm(nzInd));
        catch %#ok<CTCH>
            disp(['Method failed to determine a threshold (minmax, noise ' ...
                num2str(filterNoise) ', background ' ...
                num2str(filterBackground) '). Exiting code.'])
            return
        end
    case 'prct'
        try
            level = prctile(imageFilteredMinusBackgroundNorm(nzInd),methodValue);
        catch %#ok<CTCH>
            disp(['Method failed to determine a threshold (prct, noise ' ...
                num2str(filterNoise) ', background ' ...
                num2str(filterBackground) '). Exiting code.'])
            return
        end
    case 'user'
        try
            level = methodValue;
        catch %#ok<CTCH>
            disp(['Method failed to determine a threshold (user, noise ' ...
                num2str(filterNoise) ', background ' ...
                num2str(filterBackground) '). Exiting code.'])
            return
        end
end

%threshold the image
imageThresholded = im2bw(imageFilteredMinusBackgroundNorm,level); 

%fill holes in thresholded image to make continuous blobs
imageThresholdedFilled = imfill(imageThresholded,'holes');

% go over blobs and remove those with a size smaller that minSize
labels = bwlabel(imageThresholdedFilled);
labels = labels.*mask;
%Take the objects with and area bigger than minSize
stats = regionprops(logical(labels), 'Area'); 
idx = find([stats.Area] > minSize);

%mask with blobs
maskBlobs = ismember(labels, idx);


%% Local maxima
%get local maxima
locmaxvalues = locmax2d(imageFilteredMinusBackgroundNorm,3,[]);

% get position of local maxima
locmaxidx = find(locmaxvalues);

%get local maxima inside blobs
locmaxInBlobs = ismember(locmaxidx,find(maskBlobs));

%remove local maxima that are inside a blob
% this eliminates bright spot that are already included in maskBlobs
locmaxidx(locmaxInBlobs) = [];

%get values inside mask
locmaxInMask = ismember(locmaxidx,find(mask));

%get actual value of the local maxima
locmaxValuesInMask = locmaxvalues(locmaxidx(locmaxInMask));


%take the mean of the local maxima excluding outliers
 [locmaxMean, locmaxStd] = robustMean(locmaxValuesInMask,[],3);
 
%take values above the third std 
threshValue = locmaxMean+3*locmaxStd;

%threshold local maxima
thresh = im2bw(locmaxvalues,threshValue);

%include blobs and local maxima in the same mask
threshidx= find(thresh);
maskBlobs(threshidx) = 1;
maskBlobs = maskBlobs.*mask; 
SE = strel('disk',1);
maskBlobsVis = imdilate(maskBlobs,SE); 
 %get labels of mask with segementation and local maxima
 labels = bwlabel(maskBlobs);


%% Plotting

 if plotRes
    
    imageScaled = (image - prctile(image(:),1)) / (prctile(image(:),99) - prctile(image(:),1));
    imageScaled(imageScaled<0) = 0;
    imageScaled(imageScaled>1) = 1;
    
    %get the blob edges from the final blob mask
    SE = strel('square',3);
    edgesBlobs = imdilate(maskBlobs,SE) - maskBlobs;
    
    %give the edge pixels a value of zero in the original image
    imageScaled(edgesBlobs==1) = 0;
    
    
    %construct a 3-layered image to show blob edges on top of
    %original image
    image3Color = repmat(imageScaled,[1 1 3]);
    image3Color(:,:,1) = image3Color(:,:,1) + edgesBlobs;
    
         
    figure('Name',['segmentation_' thresholdMethod '_noise' num2str(filterNoise)...
       '_background' num2str(filterBackground) '_min=' num2str(minSize)]);
        imshow(image3Color,[]);
       
   
 

end
