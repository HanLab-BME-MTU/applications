function [ N,X,sp, rosin, otsu ] = histogramAndThresh( varargin )
%histogramAndThresh Run histogram and plot thresholds

optimalHistogram(varargin{:});

rosin = thresholdRosin(varargin{1});
otsu = thresholdOtsu(varargin{1});

line([rosin rosin],get(gca,'YLim'),'color','r');
line([otsu otsu],get(gca,'YLim'),'color','g');

legend({'Histogram','Rosin Threshold','Otsu Threshold'})

end

