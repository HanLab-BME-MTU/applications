function [stepSize,fracSteps] = getPhotoBleachStepSize(tracksFinal,numStepsRange,alpha,doPlot,checkFit,minLength)
%GETPHOTOBLEACHSTEPSIZE determines photobleaching step size by fitting multi-step functions to single-molecule tracks
%
%SYNOPSIS [stepSize,fracSteps] = getPhotoBleachStepSize(tracksFinal,numStepsRange,alpha,doPlot,checkFit,minLength)
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
%       minLength    : Minimum length of track to use for analysis.
%                      Optional. Default: 25.
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

if nargin < 6 || isempty(minLength)
    minLength = 25;
end

%% Calculation

%get amplitude matrix from tracksFinal
tracksMat = convStruct2MatIgnoreMS(tracksFinal);
ampMat = tracksMat(:,4:8:end)';
[numFrames,numTracks] = size(ampMat);

%determine which tracks end before end of movie
trackSEL = getTrackSEL(tracksMat);
trackEndsBeforeLast = trackSEL(:,2) < numFrames;

%retain only tracks longer than minLength for analysis
trackLength = zeros(numTracks,1);
for iTrack = 1 : numTracks
    trackLength(iTrack) = length(find(~isnan(ampMat(:,iTrack))));
end
indxGood = find(trackLength>=minLength);
keepTrack = zeros(numTracks,1);
keepTrack(indxGood) = 1;

numSteps = zeros(numTracks,1);
stepSize = NaN(numTracks,numStepsRange(2)+1);
for iTrack = indxGood'
    
    %fit amplitude time series
    [stepX,valY] = fitDataWithMultiStepFunc((1:numFrames)',ampMat(:,iTrack),numStepsRange,alpha,0);
    
    %get step sizes and number of steps
    numStepsTmp = length(stepX);
    numSteps(iTrack) = numStepsTmp + trackEndsBeforeLast(iTrack);
    stepSize(iTrack,max(1,numStepsTmp+trackEndsBeforeLast(iTrack))) = valY(end); %NOTE 1: tracks with zero steps get a value to avoid crashing, these will be removed later
    stepSize(iTrack,1:numStepsTmp) = abs(diff(valY))';
    
    %check fit if requested
    if checkFit
        x = (1:numFrames)';
        y = ampMat(:,iTrack);
        indxFirst = find(~isnan(y),1,'first');
        indxLast = find(~isnan(y),1,'last');
        x = x(indxFirst:indxLast);
        y = y(indxFirst:indxLast);
        if trackEndsBeforeLast(iTrack)
            stepX(end+1) = indxLast+1; %#ok<AGROW>
            valY(end+1) = 0; %#ok<AGROW>
        end
        overlayMultiStepFunc(x,y,stepX,valY,['Track ' num2str(iTrack)])
        userEntry = input(['Step Detection for Track ' num2str(iTrack) ' OK? y/n '],'s');
        close
        if strcmp(userEntry,'n')
            keepTrack(iTrack) = 0;
        end
    end
    
end

%keep info for good tracks only
indxGood = find(keepTrack);
numSteps = numSteps(indxGood);
stepSize = stepSize(indxGood,:);

%for tracks with zero steps, remove step value (see NOTE 1 above)
stepSize(numSteps==0,1) = NaN;

%re-format for output
stepSize = stepSize(~isnan(stepSize));
fracSteps = hist(numSteps,numStepsRange(1):numStepsRange(2));
fracSteps = fracSteps / sum(fracSteps);
fracSteps = [(numStepsRange(1):numStepsRange(2))' fracSteps'];

%% Plotting

if doPlot
    
    figure
    subplot(1,2,1)
    optimalHistogram(stepSize,[],0)
    title('Histogram of detected step sizes')
    subplot(1,2,2)
    plot(fracSteps(:,1),fracSteps(:,2))
    title('Fraction of tracks with specified number of steps')
    
end


%% ~~~ the end ~~~

