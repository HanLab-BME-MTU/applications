function [xPlot, yPlot, sigma, est] = ptEstimateLine (xdata, ydata, subPlot, drugTimepoint)
% ptEstimateLine plots the estimated function calculated by ptFitLine 
%
% SYNOPSIS       ptEstimateLine (xdata, ydata, slope, subPlot, drugTimepoint)
%
% INPUT          xdata   : x-axis values
%                ydata   : y-axis values
%                subPlot : number of subplot graph should be plotted in.
%                          Format is [a b c], where a is subplot row, b 
%                          subplot column and c the position. 
%                          [] if not needed
%                drugTimepoint : time point of drug application
%                
% OUTPUT         yPlot : y-axis values 
%                est : the estimated coefficients from ptFitCurve
%
% DEPENDENCIES   ptEstimateLine uses { ptFitLine }
%                                  
%                ptEstimateLine is used by { ptPlotChaosStats }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Mar 05          Initial version

% Initialize yPlot
yPlot = [];

yRel = ydata;
yRel(1:drugTimepoint) = NaN;
yRelIndx = find(~isnan(yRel));
yRel = yRel(yRelIndx);
xRel = xdata(yRelIndx);

% Calculate the tanh estimate that best fits on the data
[est, sigma, error] = ptFitLine(xRel, yRel);

% Check whether a curve could be fitted at all. If not, do not plot
if error   
   return;
end

% Calculate function values using the tanh function
xPlot = xRel;
yPlot = est(1) + est(2) * xRel;

% Plot values on the graph currently open
if ~isempty(subPlot)
    if length(subPlot) == 3
        subplot (subPlot(1),subPlot(2),subPlot(3)); 
        plot(xPlot, yPlot,'r-');
    else
        fprintf(1,'ptPlotEstimate: estimate can not be plotted in subplot.\n');
    end
else
    plot(xPlot, yPlot, 'r-');
end

% Calculate t0 and t1
t0 = xRel(1);
t1 = xRel(end);

% Calculate min and max values
minValue = yPlot(1);
maxValue = yPlot(end);

% Position of the text on the graph
xValue = max(xdata)-max(xdata)*0.3; yValue = max(ydata);
horAlign = 'right'; vertAlign = 'top';  % Alignment values

% Calculate slope value 
slopeValue = (yPlot(end) - yPlot(1)) / (xPlot(end) - xPlot(1));

% Generate the text
[lineText, errmsg] = sprintf('Min = %f\nMax = %f\nt0 = %d\nt1 = %d\nSlope = %f', ...
                             minValue, maxValue, t0, t1, slopeValue);

% Print the text on the graph
if isempty(errmsg)
   %text (xValue, yValue, lineText, 'HorizontalAlignment',horAlign,'VerticalAlignment',vertAlign);
end
