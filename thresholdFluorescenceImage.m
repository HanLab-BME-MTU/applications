function varargout = thresholdFluorescenceImage(imageIn,showPlots)

%
% mask = thresholdFluorescenceImage(imageIn)
% 
% [mask,thesholdValue] = thresholdFluorescenceImage(imageIn,showPlots)
% 
% This function thresholds the input fluorescence image by automatically
% selecting a threshold based on the intensity histogram.
% 
% Input:
% 
%   imageIn - The N-Dimensional image to be thresholded.
% 
% 
%   showPlots - If true, a plot of the histogram and an overlay of the mask
%   on the image will be shown. The overlay plot only works if the image is
%   2D.
% 
% 
% Output:
% 
%   mask - The logical array which is true at values in imageIn greater
%   than the selected threshold.
% 
%   thresholdValue - The intensity value selected for thresholding.
%
%
%Hunter Elliott, 11/7/08
%

if nargin < 2 || isempty(showPlots)
    showPlots = 0;
end

%Convert to double if necessary
imageIn = double(imageIn);

%Get histogram, using Jonas' automatic bin number selection & smoothing
[vals,bins,histSpline] = histogram(imageIn(:),'smooth');

%Find the location of extrema in the histogram
histExtrema = fnzeros(fnder(histSpline));

%Get rid of the 'fake' extrema sometimes produced at beginning and end of
%distrubution by spline form.
histExtrema = histExtrema(1,:); %remove the intervals
histExtrema = histExtrema((histExtrema ~= ... %These will always be at first or last breaks in spline
            histSpline.breaks(1)) & (histExtrema ~= histSpline.breaks(end)));
histExtVals = fnval(histSpline,histExtrema);


%Determine whether each extrema is maximum or minimum
isMax = fnval(fnder(histSpline,2),histExtrema) < 0;


%Fint the lowest-intensity maximum, assume this is the background peak.
iBackMax = find(isMax,1,'first');

%Find the first minimum after this maximum. This is used as the threshold.
iSep = iBackMax + 1;

if iSep > length(histExtrema);
    error('Could not automatically determine a threshold value!');
end

thresholdValue = histExtrema(iSep);
minVal = histExtVals(iSep);

imageMask = imageIn >= thresholdValue;

if showPlots    
    histFig = figure;
    fnplt(histSpline)    
    hold on
    plot(histExtrema,histExtVals,'ok')
    plot(thresholdValue,minVal,'xr')
    
    if ndims(imageIn) == 2    
        maskFig = figure;
        imagesc(imageIn);
        hold on
        contour(imageMask,'w')
        colormap hot    
    end
end


if nargout == 2;
    varargout{1} = imageMask;
    varargout{2} = thresholdValue;
else
    varargout = {imageMask};
end