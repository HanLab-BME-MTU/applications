function [yPlot, est] = ptPlotEstimate (xdata, ydata, slope, subPlot)
% ptPlotEstimate plots the estimated function calculated by ptFitCurve 
%
% SYNOPSIS       ptPlotEstimate (xdata, ydata)
%
% INPUT          xdata   : x-axis values
%                ydata   : y-axis values
%                slope   : expected slope (<0 neg, >0 pos)
%                subPlot : number of subplot graph should be plotted in.
%                          Format is [a b c], where a is subplot row, b 
%                          subplot column and c the position. 
%                          [] if not needed
%                
% OUTPUT         yPlot : y-axis values 
%                est : the estimated coefficients from ptFitCurve
%                (plots are also shown on the screen and saved on the disk) 
%
% DEPENDENCIES   ptPlotEstimate uses { ptFitCurve }
%                                  
%                ptPlotEstimate is used by { ptPlotNeighbourTraj
%                                            ptPlotNeighbourChanges }
%
% Revision History
% Name                  Date            Comment
% --------------------- --------        --------------------------------------------------------
% Andre Kerstens        Oct 04          Initial version
% Andre Kerstens        Dec 04          Added slope parameter to be able to
%                                       handle pos and neg slopes
% Andre Kerstens        Mar 05          Changed default fit function to
%                                       tangens hyperbolicus (tanh) 

% Initialize yPlot
yPlot = [];

% Calculate the estimate that best fits on the data
[est, error] = ptFitCurve(xdata, ydata, 35, slope);

% Check whether a curve could be fitted at all. If not, do not plot
if error   
   return;
end

% if uint16(est(3)) > uint16(est(4))
%     fprintf(1,'ptPlotEstimate: Error - Estimated values cannot be used for plotting.\n');
%     return;
% end

% Calculate function values using the tanh function
xPlot = xdata;
yPlot = est(1) + est(2) * tanh(est(3)*(xPlot-est(4)));

% Calculate t0 and t1 (take t0 to be 5% less or more than the value of the
% first plateau and t1 to be 5% less or more than the value of the second
% plateau)
if slope >= 0
    crossMin = yPlot(1) + (yPlot(end) - yPlot(1))/20; % /20 makes 5%
    crossMax = yPlot(end) - (yPlot(end) - yPlot(1))/20;
    t0 = find(yPlot>crossMin,1,'first');
    t1 = find(yPlot<crossMax,1,'last');
else
    crossMin = yPlot(1) - (yPlot(1) - yPlot(end))/20;
    crossMax = yPlot(end) + (yPlot(1) - yPlot(end))/20;    
    t0 = find(yPlot<crossMin,1,'first');
    t1 = find(yPlot>crossMax,1,'last');
end

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

% Put some text giving t0, t1 and the slope (a2-a1)/(t1-t0)
% To do this calculate some values first
if slope > 0
    maxValue = est(1)-est(2);   % value of the second plateau
    minValue = 2*est(2)+(est(1)-est(2));  % value of the 1st plateau
    
    % Position of the text on the graph
    xValue = max(xdata)-max(xdata)*0.9; yValue = max(ydata);
    horAlign = 'left'; vertAlign = 'top';  % Alignment values
else
    maxValue = 2*est(2)+(est(1)-est(2));
    minValue = est(1)-est(2);
    xValue = max(xdata)-max(xdata)*0.9; yValue = max(ydata)-max(ydata)*0.9;
    horAlign = 'left'; vertAlign = 'bottom';
end

% Calculate slope value (a value 1 for est(3) is equal to pi/180, we scale
% using est(2)
slopeValue = (est(3)/(pi/180))/est(2);

% Generate the text
[lineText, errmsg] = sprintf('Min = %f\nMax = %f\nt0 = %d\nt1 = %d\nSlope = %f', ...
                             minValue, maxValue, t0, t1, slopeValue);

% Print the text on the graph
if isempty(errmsg)
   text (xValue, yValue, lineText, 'HorizontalAlignment',horAlign,'VerticalAlignment',vertAlign);
end
