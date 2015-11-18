function [dataFit] = plotMultipleSmoothingSpline(outputDir, data, times, names, colors, plotTitle, yLabelName, smoothingPara, yMax, yMin, alignTimes)
%Scatter plots given data sets in one plot and fits smoothing spline through the scatterplots
%
%SYNOPSIS [dataFit] = plotMultipleSmoothingSpline(outputDir, data, times, names, colors, plotTitle, yLabelName, smoothingPara)
%
%INPUT
%   outputDir       : directory where the figure will be saved
%   data            : cell array of data sets (each set must be a column)
%   times           : cell array of times corresponding to the dataset (each set must be a column)
%   names           : cell array of names corresponding to each set (used in the legend)
%   colors          : cell array of colors to be used for each set
%   plotTitle       : name of the plot
%   unit            : y axis name
%   smoothingPara   : smoothing parameter for smoothing spline fit
%   yMax            : y axis range maximum
%   yMax            : y axis range minimum
%   alignTimes      : time point where the data was aligned for each
%                     condition
%   
%OUTPUT
%   dataFit         : cell array of fitObjects from the smoothing spline fit
%
%Tae H Kim, July 2015

%% Initialization
%assign default value
if isempty(smoothingPara)
    smoothingPara = .05;
end
%creates figure and stores the figure handle
figureHandle = figure('Name', plotTitle);
try
    %% Smoothing data
    %remove nan and inf and inf
    nCond = numel(data);
    for iCond = 1:nCond
        mask = ~(isnan(data{iCond}) | isinf(data{iCond}));
        data{iCond} = data{iCond}(mask);
        times{iCond} = times{iCond}(mask);
    end
    %smoothing
    smoothData = cellfun(@(x) smooth(x, 5), data, 'UniformOutput', false);
    smoothData = cellfun(@(x) x(3:end-2), smoothData, 'UniformOutput', false);
    smoothTimes = cellfun(@(x) smooth(x, 5), times, 'UniformOutput', false);
    smoothTimes = cellfun(@(x) x(3:end-2), smoothTimes, 'UniformOutput', false);
    %% Fitting
    dataFit = cellfun(@(x, y) fit(x, y, 'smoothingspline', 'smoothingParam', smoothingPara), smoothTimes, smoothData, 'UniformOutput', false);
    %% Plotting
    
    hold on;
    %plots all data and stores all line handles
    lineHandle = cellfun(@plot, dataFit, times, data, 'UniformOutput', false);
    %plot vertical lines indicating aligning Times
    nCond = numel(data);
    for iCond = 1:nCond
        plot([alignTimes(iCond), alignTimes(iCond)], [yMax, yMin], 'Color', colors{iCond});
    end
    %change the color so that color of data and fit match
    cellfun(@(x, y) set(x, 'Color', y), lineHandle, colors);
    %create legends only contain the fit
    fitHandle = [lineHandle{:}];
    legend(fitHandle(2,:), names);
    %axis limit
    if ~isempty(yMax)
        ylim([yMin, yMax]);
    end
    %label axis
    xlabel('Time (min)');
    ylabel(yLabelName);
    title(plotTitle);
    %% Saving
    %save and close
    savefig(figureHandle, [outputDir filesep plotTitle '.fig']);
    close(figureHandle);
catch err
    warning(['Could not plot ' plotTitle]);
    disp(getReport(err));
    dataFit = [];
    close(figureHandle);
end
end