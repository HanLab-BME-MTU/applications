function [dataFit] = calcMultipleSmoothingSpline(data, times, inOutFlag, smoothingPara)
%Scatter plots given data sets in one plot and fits smoothing spline through the scatterplots
%
%SYNOPSIS [dataFit] = plotMultipleSmoothingSpline(outputDir, data, times, names, colors, plotTitle, yLabelName, smoothingPara)
%
%INPUT
%   data            : cell array of data sets (each set must be a column)
%   times           : cell array of times corresponding to the dataset (each set must be a column)
%   plotTitle       : name of the plot
%   inOutFlag       : cell array of flags indicating inlier data points
%                     (value 1) and outlier data points (value 0)
%   smoothingPara   : smoothing parameter for smoothing spline fit, see fitoptions('smoothingspline')
%   
%OUTPUT
%   dataFit         : cell array of fitObjects from the smoothing spline fit
%
% See also csaps, fit, curvefit, smoothing-splines, splinetool, fitoptions('smoothingspline')
%Tae H Kim, July 2015
% Mark Kittisopikul, November 2015

%% Initialization
%assign default value
if isempty(smoothingPara)
    smoothingPara = .01;
end
%creates figure and stores the figure handle
    fxn = @(varargin) calcMultipleSmoothingSplinePerCondition(varargin{:}, smoothingPara);
    dataFit = cellfun(fxn, data, times, inOutFlag,'UniformOutput',false);
end
function dataFit = calcMultipleSmoothingSplinePerCondition(data, times, inOutFlag, smoothingPara)
    try
        %% Smoothing data
        %remove nan and inf
        mask = isfinite(data);
        data = data(mask);
        times = times(mask);
        inOutFlag = inOutFlag(mask);
        %KJ: discard outliers
        timesIn = times(inOutFlag);
        dataIn = data(inOutFlag);
        %smoothing
        smoothData = smooth(dataIn, 5);
        smoothData = smoothData(3:end-2);
        smoothTimes = smooth(timesIn, 5);
        smoothTimes = smoothTimes(3:end-2);

        %% Fitting
        % See csaps, fitoptions('smoothingspline')
        dataFit = fit(smoothTimes, smoothData, 'smoothingspline', ...
                'smoothingParam', smoothingPara);
    catch err
        dataFit = [];
    end
end