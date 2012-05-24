function [windowMotionType,fracMotionType,medianSpeedType,...
    medianWinMotionRange] = classifyEdgeMotion(protSamples,doPlot,...
    windowRange,frameRange)
%CLASSIFYEDGEMOTION classifies edge motion in each window based on protrusion vector samples
%
%SYNPOSIS [windowMotionType,fracMotionType,medianSpeedType,...
%    medianWinMotionRange] = classifyEdgeMotion(protSamples,doPlot,...
%    windowRange,frameRange)
%
%INPUT  protSample      : The protrusion samples as output by the windowing
%                         software.
%       doPlot          : Flag with value 1 to plot and 0 not to plot
%                         classification. Optional. Default: 0.
%       windowRange     : 2-row vector indicating range of windows (i.e.
%                         window number along cell edge) to include in
%                         analysis.
%                         Optional. Default: [1 (last window)].
%       frameRange      : 2-row vector indicating range of window frames
%                         to include in analysis.
%                         Optional. Default: [1 (last frame)].
%OUTPUT windowMotionType: 2D array same size and protSamples.avgNormal
%                         indicating motion classification, where +1 means
%                         protrusion, -1 means retraction, 0 means pause,
%                         and NaN indicates windows without a prutrusion
%                         vector.
%       fracMotionType  : Fraction of intervals classified as protrusion,
%                         retraction, pause and undetermined.
%       medianSpeedType : Median speed in each motion type (same order as
%                         fracMotionType).
%       medianWinMotionRange: Median motion range across all windows.
%
%Khuloud Jaqaman, December 2011

%% Input

if nargin < 1
    error('--classifyEdgeMotion: Please enter protrusion samples!');
end

if nargin < 2
    doPlot = 0;
end

%get number of windows and number of frames
avgNormal = protSamples.avgNormal;
[numWindows,numFrames] = size(avgNormal);

if nargin < 3 || isempty(windowRange)
    windowRange = [1 numWindows];
end

if nargin < 4 || isempty(frameRange)
    frameRange = [1 numFrames];
end

%% Classification

%for windows/frames outside of the range of interest, convert their values
%to NaN
avgNormal(1:windowRange(1)-1,:,:) = NaN;
avgNormal(windowRange(2)+1:end,:,:) = NaN;
avgNormal(:,1:frameRange(1)-1,:) = NaN;
avgNormal(:,frameRange(2):end,:) = NaN;

%assign edge position and protrusion vector standard deviations
edgePosStd = 1; %in pixels
protVecStd = sqrt(2*edgePosStd);

%initialize output
windowMotionType = NaN(size(avgNormal));
motionRange = NaN(numWindows,1);

%go over all windows
for iWin = 1 : numWindows
    
    %initialize motion classification vector
    classCurr = NaN(1,numFrames-1);
    
    %get current protrusion normal
    protNorm = avgNormal(iWin,1:end-1)';
    
    if any(~isnan(protNorm))
        
        %determine sign of each protrusion normal
        protNormSign = sign(protNorm);
        
        %extract points of switching
        protNormSignDiff = diff(protNormSign);
        switchPoints = find(protNormSignDiff~=0) + 1;
        switchPoints = [1; switchPoints; numFrames]; %#ok<AGROW>
        
        %calculate length of intervals between switching
        intervalLength = diff(switchPoints);
        
        %get number of intervals
        numIntervals = length(intervalLength);
        
        %calculate position change (i.e. cumulative protrusion/retraction) in
        %each interval
        intervalProt = NaN(numIntervals,1);
        for iInt = 1 : numIntervals
            intervalProt(iInt) = sum(protNorm(switchPoints(iInt):switchPoints(iInt+1)-1));
        end
        
        %initialize interval classification vector
        intervalClass = NaN(size(intervalLength));
        
        %assign as significant any intervals that are >= 3 frame steps
        indx = find(intervalLength >= 3);
        intervalClass(indx) = sign(intervalProt(indx)); % +1 = protrusion, -1 = retraction
        
        %go over intervals that are <= 2 frame steps and test for significance
        indx2 = find(intervalLength<=2);
        oneMinusPvalue = normcdf(abs(intervalProt(indx2)),0,protVecStd);
        indx = indx2(oneMinusPvalue>0.95);
        intervalClass(indx) = sign(intervalProt(indx));% +1 = protrusion, -1 = retraction
        indx = indx2(oneMinusPvalue<=0.95);
        intervalClass(indx) = 0; % 0 = pause
        
        %         %assign as significant any intervals that are >= 2 frame steps
        %         indx = find(intervalLength >= 2);
        %         intervalClass(indx) = sign(intervalProt(indx)); % +1 = protrusion, -1 = retraction
        %
        %         %go over intervals that are 1 frame step and test for significance
        %         indx2 = find(intervalLength==1);
        %         oneMinusPvalue = normcdf(abs(intervalProt(indx2)),0,protVecStd);
        %         indx = indx2(oneMinusPvalue>0.95);
        %         intervalClass(indx) = sign(intervalProt(indx));% +1 = protrusion, -1 = retraction
        %         indx = indx2(oneMinusPvalue<=0.95);
        %         intervalClass(indx) = 0; % 0 = pause
        
        %from intervalClass, assign classification to individual frame steps
        for iInt = 1 : numIntervals
            classCurr(switchPoints(iInt):switchPoints(iInt+1)-1) = intervalClass(iInt);
        end
        
        %save in output variable
        windowMotionType(iWin,1:end-1) = classCurr;
        
        %integrate protrusion normals to get position over time for
        %plotting
        posNorm = nancumsumwithzero(protNorm);
        
        %calculate edge range of motion
        motionRange(iWin) = max(posNorm)-min(posNorm);
        
        if doPlot
            
            %plot protrusion in green, retraction in red and pause in cyan
            %intervals with missing points stay empty
            figure, hold on
            plot(posNorm,'kd'), myErrorbar(posNorm,edgePosStd*ones(numFrames,1))
            indx = find(classCurr==1);
            for iInt = indx
                plot([iInt iInt+1],posNorm(iInt:iInt+1),'g')
            end
            indx = find(classCurr==-1);
            for iInt = indx
                plot([iInt iInt+1],posNorm(iInt:iInt+1),'r')
            end
            indx = find(classCurr==0);
            for iInt = indx
                plot([iInt iInt+1],posNorm(iInt:iInt+1),'c')
            end
            
        end %(if doPlot)
        
    end %(if any(~isnan(protNorm)))
    
end %(for iWin = 1 : numWindows)

%% Edge activity characteristics

%find windows in each category
indxProt = find(windowMotionType==1);
indxRet = find(windowMotionType==-1);
indxPause = find(windowMotionType==0);
indxUndet = find(isnan(windowMotionType));

%get fraction of time spent in each category
fracMotionType = [length(indxProt) length(indxRet) length(indxPause) length(indxUndet)] / length(windowMotionType(:));

%get median speed in each category
medianSpeedType = [median(avgNormal(indxProt)) median(avgNormal(indxRet)) median(avgNormal(indxPause)) median(avgNormal(indxUndet))];

%get median motion range for all windows
medianWinMotionRange = nanmedian(motionRange);

%% ~~~ the end ~~~
