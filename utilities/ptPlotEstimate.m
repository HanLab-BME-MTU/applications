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

% Initialize yPlot
yPlot = [];

% Calculate the estimate that best fits on the data
[est, error] = ptFitCurve(xdata,ydata,35,slope);

% Check whether a curve could be fitted at all. If not, do not plot
if error   
   return;
end

if uint16(est(3)) > uint16(est(4))
    fprintf(1,'ptPlotEstimate: Error - Estimated values cannot be used for plotting.\n');
    return;
end

% Chop x-axis data up into the 3 estimated parts
xLeft = xdata(1:int16(est(3))-1);
xMiddle = xdata(int16(est(3)):int16(est(4))-1);
xRight = xdata(int16(est(4)):end);

% Put all the x values together again
xPlot = [xLeft xMiddle xRight];

% Chop y-axis data up into the 3 estimated lines
yLeft = ones(1,length(xLeft))*est(1);
yMiddle = ((est(2)-est(1))/(est(4)-est(3)))*(xMiddle-est(3)) + est(1);
yRight = ones(1,length(xRight))*est(2);

% Put the y values together
yPlot = [yLeft yMiddle yRight];

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
[lineText, errmsg] = sprintf('t0 = %d\nt1 = %d\nslope = %f', int16(est(3)), ...
                             int16(est(4)), ((est(2)-est(1))/(est(4)-est(3))));
if slope > 0
    xValue = max(xdata)-max(xdata)*0.9; yValue = max(ydata);
    horAlign = 'left'; vertAlign = 'top';
else
    xValue = max(xdata)-max(xdata)*0.9; yValue = max(ydata)-max(ydata)*0.9;
    horAlign = 'left'; vertAlign = 'bottom';
end

text (xValue, yValue, lineText, 'HorizontalAlignment',horAlign,'VerticalAlignment',vertAlign);
