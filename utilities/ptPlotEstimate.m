function [yPlot, est] = ptPlotEstimate (xdata, ydata, slope, subPlot, drugTimepoint)
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
%                drugTimepoint : time point of drug application
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

% Calculate the tanh estimate that best fits on the data
[est, dataUsed, sigma, error] = ptFitCurve(xdata, ydata, drugTimepoint, slope);

% Check whether a curve could be fitted at all. If not, do not plot
if error   
   return;
end

% Calculate function values using the tanh function
xPlot = xdata;
yPlot = est(1) + est(2) * tanh(est(3)*(xPlot-est(4)));

% Determine first estimate of t0 and t1
t0 = round(est(4) - abs(1/est(3)));
t1 = round(est(4) + abs(1/est(3)));

% Use these estimates to find find values based on the curve fit error
% (which basically means we take into account the amount of variation of
% the data itself)
% First determine whether there is a plateau at the beginning of the curve
if t0 <= 1 | t0 >= length(xPlot) % No plateau present
    t0 = 0;
else  % A plateau exists
    % extract the relevant part of the data
    yRel = yPlot(1:t0);
    dataUsedRel = dataUsed(find(dataUsed<=t0));
    yRel = yRel(dataUsedRel);
        
    % Calculate average error between data and fitted values
    error = sqrt((ydata(dataUsedRel) - yRel).^2);
    avgError = sum(error)/length(error);
    
    % Determine the point where the fitted curve has a distance >= avgError
    % This is our new t0
    if est(3) > 0
        t0 = find(yPlot >= (yPlot(1)+avgError),1,'first');
    else
        t0 = find(yPlot <= (yPlot(1)-avgError),1,'first');
    end
    
    % There are cases where t0 is [] because the error is so large. In that
    % case stick to the original estimate
    if isempty(t0)
        t0 = round(est(4) - abs(1/est(3)));
    end
end

% Do the same for the plateau at the end
if t1 >= length(xPlot) | t1 <= 1 % No plateau present
    t1 = 0;
else  % A plateau exists
    % extract the relevant part of the data
    yRel = yPlot(round(t1):end);
    dataUsedRel = dataUsed(find(dataUsed>t1));
    dataUsedRelEnd = dataUsedRel - t1;
    yRel = yRel(dataUsedRelEnd);
        
    % Calculate average error between data and fitted values
    error = sqrt((ydata(dataUsedRel) - yRel).^2);
    avgError = sum(error)/length(error);
    
    % Determine the point where the fitted curve has a distance <= avgError
    % This is our new t1
    if est(3) > 0
        t1 = find(yPlot <= (yPlot(end)-avgError),1,'last');
    else
        t1 = find(yPlot >= (yPlot(end)+avgError),1,'last');
    end
    
    % There are cases where t0 is [] because the error is so large. In that
    % case stick to the original estimate. Same for cases where t1 is
    % smaller than t0
    if (isempty(t1)) | (t1 < t0)
        t1 = round(est(4) + abs(1/est(3)));
    end
end
   
% Last check to make sure t1 is larger than t0
if t1 < t0
    t1 = t0;
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

% Calculate min and max values and alignment values of the text
if est(3) > 0        
    % Calculate min and max values
    minValue = yPlot(1);
    maxValue = yPlot(end);
    
    % Position of the text on the graph
    xValue = max(xdata)-max(xdata)*0.9; yValue = max(ydata);
    horAlign = 'left'; vertAlign = 'top';  % Alignment values
else    
    % Calculate min and max values
    minValue = yPlot(end);
    maxValue = yPlot(1);
    
    xValue = max(xdata)-max(xdata)*0.1; yValue = max(ydata);
    horAlign = 'right'; vertAlign = 'top';
end

% Calculate slope value (a value 1 for est(3) is equal to pi/180, we scale
% using est(2)
%slopeValue = (est(3)/(pi/180))/est(2);
if t0 ~= 0 & t1 ~= 0
    if (t1 - t0) ~= 0
        slopeValue = (yPlot(t1) - yPlot(t0)) / (t1 - t0);
    else
        slopeValue = (yPlot(t1) - yPlot(t0));
    end
else
    slopeValue = 0;
end

% Formatting related stuff
if t0 == 0
    t0 = 'N/A';
    t0Format = '%s';
else
    t0Format = '%d';
end
if t1 == 0
    t1 = 'N/A';
    t1Format = '%s';
else
    t1Format = '%d'; 
end
if slopeValue == 0
    slopeValue = 'N/A';
    slopeFormat = '%s';
else
    slopeFormat = '%f';
end

% Generate the text
[lineText, errmsg] = sprintf(['Min = %f\nMax = %f\nt0 = ' t0Format '\nt1 = ' ...
                              t1Format '\nSlope = ' slopeFormat '\nSigma = %f\n'], ...
                              minValue, maxValue, t0, t1, slopeValue, sigma);    

% Print the text on the graph
if isempty(errmsg)
   text (xValue, yValue, lineText, 'HorizontalAlignment',horAlign,'VerticalAlignment',vertAlign);
end
