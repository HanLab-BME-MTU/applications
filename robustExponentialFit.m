function [u,sigmaU] = robustExponentialFit(xData, yData, fitConstant, verbose)
%ROBUSTEXPONENTIALFIT uses least median squares to fit an exponential model
%
% SYNOPSIS [u,sigmaU] = robustExponentialFit(xData, yData, fitConstant)
%
% INPUT    xData : x-values of the data.
%          yData : y-values of the data.
%             xData and yData have to be vectors of the same length
%          fitConstant (opt): default: 0
%               - 0: model will be y = A*exp(-x/B)
%               - 1: model will be y = A*exp(-x/B) + C
%          verbose (opt): Whether to plot a figure with the results [{0}/1]
%
% OUTPUT   u: vector of the fitted values A,B and potentially C
%          sigmaU: a-posteriori uncertainty of the parameters
%
% c: jonas 11/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%================
% TEST INPUT
%================

% defaults
def_fitCon = 0;
def_verbose = 0;

% input arguments
if nargin < 2 || isempty(xData) || isempty(yData)
    error('please specify at least 2 input arguments!')
end
if length(xData) ~= length(yData)
    error('xData and yData need to have the same length!')
end

% make sure the vectors are col-vectors
xData = returnRightVector(xData);
yData = returnRightVector(yData);

% fitConstant
if nargin < 3 || isempty(fitConstant)
    fitConstant = def_fitCon;
elseif numel(fitConstant)>1 || ~any(ismember(fitConstant,[0,1]))
    error('fitConstant should either be 0 or 1!')
end

% verbose
if nargin < 4 || isempty(verbose)
    verbose = def_verbose;
elseif numel(verbose)>1 || ~any(ismember(verbose,[0,1]))
    error('verbose should either be 0 or 1!')
end

%===================


%==================
% EXPONENTIAL FIT
%==================

% get initial guess
init = zeros(1,3);
if fitConstant
    init(3) = robustMean(yData(end-4:end)); %C
else
    init(3) = 0;
end
init(1) = robustMean(yData(1:5))-init(3); %A
% get init2 from linear fit
ylin = log((yData-init(3))./init(1));
% real prevents complex numbers (we only initialize here!)
init(2) = -robustMean(xData./real(ylin));


parameters = struct('xdata',xData,'ydata',yData);
options = optimset('Display','off');

if fitConstant
    [u,dummy,sigmaU]=leastMedianSquares('(ydata-(u(1)*exp(-xdata/u(2)))-u(3))',init,options,parameters);
else
    [u,dummy,sigmaU]=leastMedianSquares('(ydata-(u(1)*exp(-xdata/u(2))))',init(1:2),options,parameters);
end

% plot results
if verbose
    fh = figure;
    % plot fit to data
    subplot(2,1,1);
    yinit = init(1)*exp(-xData./init(2))+init(3);
    plot(xData,yData,'.', xData,yinit,'b');

    figure(fh)
    hold on
    if fitConstant
        yfit = u(1)*exp(-xData./u(2))+u(3);
        plot(xData,yfit,'r')
    else
        yfit = u(1)*exp(-xData./u(2));
        plot(xData,yfit,'r')
    end
    % plot residuals
    subplot(2,1,2)
    stem(xData, yData-yinit,'.b')
    hold on
    stem(xData,yData-yfit,'.r')
end
