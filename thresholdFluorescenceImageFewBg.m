function thresholdValue = thresholdFluorescenceImageFewBg(imageIn,showPlots)

%
% thresholdValue = thresholdFluorescenceImage(imageIn)
% 
% thresholdValue = thresholdFluorescenceImage(imageIn,showPlots)
% 
% This function selects a threshold for the input fluorescence image by
% analyzing the image's intensity distribution. This requires good signal-
% to-noise. It only works if there is few BG in the image.
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
% Output:
% 
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
[~,~,histSpline] = histogram(imageIn(:),'smooth');

%Find the location of extrema in the histogram
histExtrema = fnzeros(fnder(histSpline));

histExtrema = histExtrema(1,:); %remove the intervals
histExtVals = fnval(histSpline,histExtrema);

%Determine whether each extrema is maximum or minimum
% isMax = fnval(fnder(histSpline,2),histExtrema) < 0;
% isMin = ~isMax;

% get the lowest minimum in the lowest quartile of the data:
cutOff = quantile(imageIn(:),0.25);
checkVec=(histExtrema<cutOff);
[~,iSep]=min(histExtVals(checkVec));

thresholdValue=histExtrema(iSep);
minVal = histExtVals(iSep);



% %Find the lowest-intensity maximum, assume this is the background peak.
% iBackMax = find(isMax,1,'first');
% 
% %Find the first minimum after this maximum. This is used as the threshold.
% iSep = iBackMax + 1;
% 
% % we don't pick minima where the difference between background maximum and
% % minimum is marginal:
% if nargin>2 && noisy && numel(histExtVals)>20
%     extVal= fnval(histSpline,histExtrema);
%     extDiff=(extVal(2:end)-extVal(1:end-1));
%     madDiff=mad(abs(extDiff));
%     
%     while iBackMax<=length(histExtrema) && (extVal(iBackMax)-extVal(iBackMax+1))<madDiff
%         iBackMax=iBackMax+2;
%         iSep = iBackMax + 1;
%     end
% end
% 
% if iSep > length(histExtrema);
%     error('Could not automatically determine a threshold value!');
% end
% 
% thresholdValue = histExtrema(iSep);
% minVal = histExtVals(iSep);

if showPlots    
    imageMask = imageIn >= thresholdValue;
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
