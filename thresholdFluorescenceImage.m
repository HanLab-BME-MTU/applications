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
% NOTE: I'm not really sure why this function works. In some rare cases it
% works as intended. Most of the time, a combination of errors
% adds up to it somehow working anyways...
%
%Hunter Elliott, 11/7/08
%

if nargin < 2 || isempty(showPlots)
    showPlots = 0;
end


%Convert to double if necessary
imageIn = cast(imageIn,'double');

%Get the distribution of values in spline form
[vals,bins,histSpline] = histogram(imageIn,'smooth');


%take the derivative of the distribution
histDeriv = fnder(histSpline);

%Find the location of extrema in the histogram
histExtrema = fnzeros(histDeriv);

%Get rid of the 'fake' extrema sometimes produced at beginning and end of
%distrubution by spline form.

minFrac = 1e-3; %Small number for eliminating end extrema
histExtrema = histExtrema(1,:); %Dump the intervals
histExtrema = histExtrema(histExtrema > (range(bins) * minFrac)); %Too small
histExtrema = histExtrema(histExtrema < (range(bins) * (1-minFrac))); %Too big
histExtVals = fnval(histSpline,histExtrema); %Just right... evaluate at these extrema

%Check number of remaingin extrema.
nExtrema = size(histExtrema,2);
if nExtrema < 2    
    error('Could not automatically determine a threshold value!');        
end





%Assuming the lowest-intensity maxima is the background, the first minima
%after this peak should be a good place to seperate background &
%foreground.
thresholdValue = histExtrema(2);
minVal = histExtVals(2);


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