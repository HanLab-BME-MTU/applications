function varargout = thresholdFluorescenceImage(imageIn,showPlots)

%Masks the input image by automatically choosing a threshold value. Beta
%version...  ;)

%Hunter Elliott, 11/7/08

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

nExtrema = size(histExtrema,2);

if nExtrema < 4    
    error('Could not automatically determine a threshold value!');        
end

%CHeck the intevals returned by fnzeros????!!!?? TEMP


%Get rid of the fake extrema produced at beginning and end of distrubution
%by spline form.
histExtrema = histExtrema(:,2:end-1);
histExtVals = fnval(histSpline,histExtrema(1,:));


thresholdValue = histExtrema(1,2);
minVal = histExtVals(2);


imageMask = imageIn >= thresholdValue;

if showPlots    
    histFig = figure;
    fnplt(histSpline)    
    hold on
    plot(histExtrema(1,:),histExtVals,'ok')
    plot(thresholdValue,minVal,'xr')
    
    maskFig = figure;
    imagesc(imageIn);
    hold on
    contour(imageMask,'w')
    colormap hot    
end


if nargout == 2;
    varargout{1} = imageMask;
    varargout{2} = thresholdValue;
else
    varargout = {imageMask};
end