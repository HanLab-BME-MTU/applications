function [maskComb,imageMinusBackground] = detectThreshLocMax(image,...
    thresholdMethod,methodValue,filterNoise,filterBackground,minSize,...
    alphaLocMax,plotRes,mask)
%detectThreshLocMax combines blob segmentation with local maxima detection
%
%SYNOPSIS maskComb = detectThreshLocMax(image,thresholdMethod,methodValue,...
%    filterNoise,filterBackground,minSize,alphaLocMax,plotRes,mask)
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
%       alphaLocMax: Alpha-value for judging the significance of local
%                   maxima.
%       plotRes   : 1 to plot segmentation results, 0 otherwise.
%                   Optional. Default: 0.
%       mask      : Binary mask. Optional. If not provided, the whole image
%                   domain is segmented.
%
%OUTPUT maskBlobs : Mask of blobs. 1 inside blobs, 0 outside.
%
%Khuloud Jaqaman January 2013

%% Output
maskComb = [];
imageMinusBackground = [];

%% Input

%check number of input arguments
if nargin < 1
    disp('Please enter at least image to be segmented');
    return
end

%get image size
[imageSizeX,imageSizeY] = size(image);

%thresholding method
if nargin < 2 || isempty(thresholdMethod)
    thresholdMethod = 'Otsu';
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
    mask = ones(imageSizeX,imageSizeY);
end

if ~logical(mask)
    error('Mask must be a logical image.');
end

%% Image preparation

%make sure that image is in double format
image = double(image);

%estimate background by filtering image with a wide Gaussian
if filterBackground > 0
    imageBackground = filterGauss2D(image,filterBackground);
else
    imageBackground = zeros(imageSizeX,imageSizeY);
end

%remove background from image
imageMinusBackground = image - imageBackground;

%remove noise by filtering image-background with a narrow Gaussian
if filterNoise > 0
    imageMinusBackgroundFiltered = filterGauss2D(imageMinusBackground,filterNoise);
else
    imageMinusBackgroundFiltered = imageMinusBackground;
end

%estimate noise per pixel
imageMinusBackgroundNoise = imageMinusBackground - imageMinusBackgroundFiltered;

%crop image
imageMinusBackgroundFiltered = imageMinusBackgroundFiltered .* mask;

%find nonzero values (due to masking)
nzInd = find(imageMinusBackgroundFiltered);

%get minumum and maximum pixel values in image
minSignal = min(imageMinusBackgroundFiltered(nzInd));
maxSignal = max(imageMinusBackgroundFiltered(nzInd));

%normalize nonzero value between 0 and 1
imageMinusBackgroundFilteredNorm = zeros(size(imageMinusBackgroundFiltered));
imageMinusBackgroundFilteredNorm(nzInd) = (imageMinusBackgroundFiltered(nzInd) - minSignal) / (maxSignal - minSignal);

%% Thresholding

%estimate the intensity level to use for thresholding the image
switch thresholdMethod
    case 'otsu'
        try
            level = graythresh(imageMinusBackgroundFilteredNorm(nzInd));
        catch %#ok<CTCH>
            disp(['Method failed to determine a threshold (otsu, noise ' ...
                num2str(filterNoise) ', background ' ...
                num2str(filterBackground) '). Exiting code.'])
            return
        end
    case 'rosin'
        try
            [dummy,level] = cutFirstHistMode(imageMinusBackgroundFilteredNorm(nzInd),0); %#ok<ASGLU>
        catch %#ok<CTCH>
            disp(['Method failed to determine a threshold (rosin, noise ' ...
                num2str(filterNoise) ', background ' ...
                num2str(filterBackground) '). Exiting code.'])
            return
        end
    case 'minmax'
        try
            level = thresholdFluorescenceImage(imageMinusBackgroundFilteredNorm(nzInd));
        catch %#ok<CTCH>
            disp(['Method failed to determine a threshold (minmax, noise ' ...
                num2str(filterNoise) ', background ' ...
                num2str(filterBackground) '). Exiting code.'])
            return
        end
    case 'prct'
        try
            level = prctile(imageMinusBackgroundFilteredNorm(nzInd),methodValue);
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
imageThresholded = im2bw(imageMinusBackgroundFilteredNorm,level);

%fill holes in thresholded image to make continuous blobs
imageThresholdedFilled = imfill(imageThresholded,'holes');

%go over blobs and remove those with a size smaller that minSize
labels = bwlabel(imageThresholdedFilled);
stats = regionprops(labels, 'Area'); %#ok<MRPBW>
idx = find([stats.Area] > minSize);

%output final blob mask from thresholding
maskBlobs = ismember(labels, idx);

%% Local maxima

% %estimate local background and noise statistics
% [bgMean,bgStd] = ...
%     spatialMovAveBG(imageFilteredMinusBackground,imageSizeX,imageSizeY);

% %estimate background/noise statistics
% bgIntDistr = imageFilteredMinusBackground(~maskBlobs);
% [bgMean,bgStd] = robustMean(bgIntDistr);
% bgStd = max(bgStd,eps);

%estimate background/noise statistics
bgMean = 0;
[~,bgStd] = robustMean(imageMinusBackgroundNoise(~maskBlobs));

%call locmax2d to get local maxima in filtered image
fImg = locmax2d(imageMinusBackgroundFiltered,[3 3],1);

%get positions and amplitudes of local maxima
localMax1DIndx = find(fImg);
[localMaxPosX,localMaxPosY,localMaxAmp] = find(fImg);

% %get background values corresponding to local maxima
% bgMeanMax = bgMean(localMax1DIndx);
% bgStdMax = bgStd(localMax1DIndx);

%calculate the p-value corresponding to the local maxima's amplitudes
%assume that background intensity is normally
%distributed with mean bgMeanMax and standard deviation bgStdMax
% pValue = 1 - normcdf(localMaxAmp,bgMeanMax,bgStdMax);
pValue = 1 - normcdf(localMaxAmp,bgMean,bgStd);

%calculate the threshold to distinguish significant local maxima
[~,threshLocMax] = cutFirstHistMode(localMaxAmp,0);

%retain only those maxima with significant amplitude
% indxKeep = pValue < alphaLocMax;
indxKeep = localMaxAmp > threshLocMax;
localMax1DIndx = localMax1DIndx(indxKeep);
localMaxPosX = localMaxPosX(indxKeep);
localMaxPosY = localMaxPosY(indxKeep);

%make a mask image from the local maxima
maskLocMax = zeros(imageSizeX,imageSizeY);
maskLocMax(localMax1DIndx) = 1;
SE = strel('square',3);
maskLocMax = imdilate(maskLocMax,SE);

%% Final mask

maskComb = maskBlobs | maskLocMax;

%% Plotting

if plotRes
    
    %figure 1: the different analysis steps
    figure

    %subplot 1: original image
    subplot(2,2,1)
    imshow(image,[prctile(image(:),1) prctile(image(:),99)]);
    
    %subplot 2: background image
    subplot(2,2,2)
    imshow(imageBackground,[prctile(imageBackground(:),1) prctile(imageBackground(:),99)]);
    
    %subplot 3: bandpass-filtered image
    subplot(2,2,3)
    imshow(imageMinusBackgroundFiltered,[prctile(imageMinusBackgroundFiltered(:),1) prctile(imageMinusBackgroundFiltered(:),99)])
    
    %subplot 4: blob edges + local maxima
    subplot(2,2,4)
    imshow(image,[prctile(image(:),1) prctile(image(:),99)]);
    hold on
    maskBounds = bwboundaries(maskBlobs);
    cellfun(@(x)(plot(x(:,2),x(:,1),'r','LineWidth',1)),maskBounds);
    plot(localMaxPosY,localMaxPosX,'go')
    
    %figure 2: final mask
    figure
    imshow(image,[prctile(image(:),1) prctile(image(:),99)]);
    hold on
    maskBounds = bwboundaries(maskComb);
    cellfun(@(x)(plot(x(:,2),x(:,1),'r','LineWidth',1)),maskBounds);
    
end

%% ~~~ the end ~~~

