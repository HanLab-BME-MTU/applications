function [stepSize,fracSteps] = getPhotoBleachStepSize(tracksFinal,numStepsRange,alpha,doPlot,checkFit)
%GETPHOTOBLEACHSTEPSIZE determines photobleaching step size by fitting multi-step functions to single-molecule tracks
%
%SYNOPSIS [stepSize,fracSteps] = getPhotoBleachStepSize(tracksFinal,numStepsRange,alpha,doPlot)
%
%INPUT  tracksFinal  : Tracker output.
%       numStepsRange: Row vector indicating minimum and maximum number
%                      steps to attempt to fit.
%                      Optional. Default: [0 5].
%       alpha        : Alpha-value for F-test to determine number of steps.
%                      Optional. Default: 0.05.
%       doPlot       : 1 to plot step size histogram, 0 otherwise.
%                      Optional. Default: 0.
%       checkFit     : 1 to visually check fit, 0 otherwise.
%                      Optional. Default: 0.
%
%OUTPUT stepSize     : Distribution of obtained step sizes.
%       fracSteps    : Fraction of tracks with 0 steps, 1 step, 2 steps, etc.
%
%Khuloud Jaqaman, August 2014

%% Output
stepSize = [];
fracSteps = [];

%% Input

%check whether correct number of input arguments was used
if nargin < 1
    disp('--getPhotoBleachStepSize: Incorrect number of input arguments!');
    return
end

if nargin < 2 || isempty(numStepsRange)
    numStepsRange = [0 5];
end

if nargin < 3 || isempty(alpha)
    alpha = 0.05;
end

if nargin < 4 || isempty(doPlot)
    doPlot = 0;
end

if nargin < 5 || isempty(checkFit)
    checkFit = 0;
end

%% Calculation

%get amplitude matrix from tracksFinal
tracksMat = convStruct2MatIgnoreMS(tracksFinal);
ampMat = tracksMat(:,4:8:end)';

[numFrames,numTracks] = size(ampMat);
numSteps = zeros(numTracks,1);
stepSize = NaN(numTracks,numStepsRange(2));
keepTrack = ones(numTracks,1);
for iTrack = 1 : numTracks
    
    %find number of datapoints in track; analyze only if it has more than
    %10 data points
    trackLength = length(find(~isnan(ampMat(:,iTrack))));
    if trackLength < 10
        keepTrack(iTrack) = 0;
    else
        
        %fit amplitude time series
        [stepX,valY] = fitDataWithMultiStepFunc((1:numFrames)',ampMat(:,iTrack),numStepsRange,alpha,0);
        
        numStepsTmp = length(stepX);
        if numStepsTmp > 0
            
            %get step size and number of steps
            numSteps(iTrack) = numStepsTmp;
            stepSize(iTrack,1:numStepsTmp) = abs(diff(valY))';
            
            %check fit if requested
            if checkFit
                x = (1:numFrames)';
                y = ampMat(:,iTrack);
                indxFirst = find(~isnan(y),1,'first');
                indxLast = find(~isnan(y),1,'last');
                x = x(indxFirst:indxLast);
                y = y(indxFirst:indxLast);
                overlayMultiStepFunc(x,y,stepX,valY,['Track ' num2str(iTrack)])
                userEntry = input(['Step Detection for Track ' num2str(iTrack) ' OK? y/n '],'s');
                close
                if strcmp(userEntry,'n')
                    keepTrack(iTrack) = 0;
                end
            end
            
        end
        
    end
            
end

%keep info for good tracks only
indxGood = find(keepTrack);
numSteps = numSteps(indxGood);
stepSize = stepSize(indxGood,:);

%re-format for output
stepSize = stepSize(~isnan(stepSize));
fracSteps = hist(numSteps,numStepsRange(1):numStepsRange(2));
fracSteps = fracSteps / sum(fracSteps);
fracSteps = [(numStepsRange(1):numStepsRange(2))' fracSteps'];

%% Plotting

if doPlot
    
    figure
    subplot(1,2,1)
    histogram(stepSize,[],0)
    title('Histogram of detected step sizes')
    subplot(1,2,2)
    plot(fracSteps(:,1),fracSteps(:,2))
    title('Fraction of tracks with specified number of steps')
    
end


%% ~~~ the end ~~~

