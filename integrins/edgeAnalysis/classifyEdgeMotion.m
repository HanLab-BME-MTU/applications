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
%OUTPUT windowMotionType: Size(protSamples.avgNormal)-by-6 array
%                         indicating motion classification, where 1 means
%                         protrusion, 2 means retraction, 3 means pause,
%                         and NaN indicates windows without a prutrusion
%                         vector.
%                         In the 3rd dimension:
%                         Entry 1: classification of each interval.
%                         Entry 2: classification of previous period.
%                         Note: period = group of consecutive intervals
%                         with same classification.
%                         Entry 3: classification of next period.
%                         Entry 4: duration of the period that the interal
%                         belond to.
%                         Entry 5: duration of previous period.
%                         Entry 6: duration of next period.
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
windowMotionType = NaN([size(avgNormal) 6]);
motionRange = NaN(numWindows,1);

%go over all windows
for iWin = 1 : numWindows
    
    %initialize motion classification vectors
    [classCurr,classBef,classAft,durCurr,durBef,durAft] = deal(NaN(1,numFrames-1));
    
    %get current protrusion normal
    protNorm = avgNormal(iWin,1:end-1)';
    
    if any(~isnan(protNorm))
        
        %classify each interval as protrusion, retraction or tentative pause
        %1 = protrusion
        %2 = retraction
        %3 = tentative pause
        %NaN = missing data point
        classCurr( protNorm>protVecStd ) = 1;
        classCurr( protNorm<-protVecStd ) = 2;
        classCurr( protNorm>=-protVecStd & protNorm<=protVecStd) = 3;
        
        %if a tentative pause interval has a protrusion normal > relaxFactor*protVecStd,
        %and at least one of its neighbors has a protrusion normal > relaxFactor*protVecStd,
        %convert it to a protrusion
        %do the equivalent for retraction tentative pauses
        relaxFactor = 0.5;
        indxPause = find(classCurr==3);
        indxPauseS = indxPause(indxPause==1);
        indxPauseE = indxPause(indxPause==(numFrames-1));
        indxPause = setdiff(indxPause,[indxPauseS indxPauseE]);
        indxProt = [ indxPause( protNorm(indxPause)'>relaxFactor*protVecStd & ...
            (protNorm(indxPause-1)'>relaxFactor*protVecStd | ...
            protNorm(indxPause+1)'>relaxFactor*protVecStd) ) ...
            indxPauseS( protNorm(indxPauseS)>relaxFactor*protVecStd & ...
            protNorm(indxPauseS+1)>relaxFactor*protVecStd ) ...
            indxPauseE( protNorm(indxPauseE)>relaxFactor*protVecStd & ...
            protNorm(indxPauseE-1)>relaxFactor*protVecStd ) ];
        indxRetr = [ indxPause( protNorm(indxPause)'<-relaxFactor*protVecStd & ...
            (protNorm(indxPause-1)'<-relaxFactor*protVecStd | ...
            protNorm(indxPause+1)'<-relaxFactor*protVecStd) ) ...
            indxPauseS( protNorm(indxPauseS)<-relaxFactor*protVecStd & ...
            protNorm(indxPauseS+1)<-relaxFactor*protVecStd ) ...
            indxPauseE( protNorm(indxPauseE)<-relaxFactor*protVecStd & ...
            protNorm(indxPauseE-1)<-relaxFactor*protVecStd ) ];
        classCurr(indxProt) = 1;
        classCurr(indxRetr) = 2;
        
        %determine the start and duration of each
        %protrusion/retraction/pause/undetermined period
        periodStart = find(diff([NaN classCurr])~=0);
        periodEnd   = find(diff([classCurr NaN])~=0);
        periodClass = classCurr(periodStart);
        periodDuration = periodEnd - periodStart + 1;
        numPeriod = length(periodClass);
        
        %sub-divide each classification based on what's before it and
        %what's after it, also based on its duration and the duration of
        %what's before it and what's after it
        if numPeriod == 1
            durCurr(periodStart(1):periodEnd(1)) = periodDuration(1);
        else
            classAft(periodStart(1):periodEnd(1)) = periodClass(2);
            durCurr(periodStart(1):periodEnd(1)) = periodDuration(1);
            durAft(periodStart(1):periodEnd(1)) = periodDuration(2);
            for iPeriod = 2 : numPeriod-1
                classBef(periodStart(iPeriod):periodEnd(iPeriod)) = periodClass(iPeriod-1);
                classAft(periodStart(iPeriod):periodEnd(iPeriod)) = periodClass(iPeriod+1);
                durCurr(periodStart(iPeriod):periodEnd(iPeriod)) = periodDuration(iPeriod);
                durBef(periodStart(iPeriod):periodEnd(iPeriod)) = periodDuration(iPeriod-1);
                durAft(periodStart(iPeriod):periodEnd(iPeriod)) = periodDuration(iPeriod+1);
            end
            classBef(periodStart(end):periodEnd(end)) = periodClass(end-1);
            durCurr(periodStart(end):periodEnd(end)) = periodDuration(end);
            durBef(periodStart(end):periodEnd(end)) = periodDuration(end-1);
        end
        
        %save in output variable
        windowMotionType(iWin,1:end-1,1) = classCurr;
        windowMotionType(iWin,1:end-1,2) = classBef;
        windowMotionType(iWin,1:end-1,3) = classAft;
        windowMotionType(iWin,1:end-1,4) = durCurr;
        windowMotionType(iWin,1:end-1,5) = durBef;
        windowMotionType(iWin,1:end-1,6) = durAft;
        
        %integrate protrusion normals to get position over time for
        %plotting
        posNorm = nancumsumwithzero(protNorm);
        
        %calculate edge range of motion
        motionRange(iWin) = max(posNorm)-min(posNorm);
        
        if doPlot
            
            %show plot
            %color-coding:
            %green  = protrusion (1)
            %red    = retraction (2)
            %cyan   = pause (3)
            %intervals with missing points stay empty
            figure, hold on
            plot(posNorm,'kd'), myErrorbar(posNorm,edgePosStd*ones(numFrames,1))
            indx = find(classCurr==1);
            for iInt = indx
                plot([iInt iInt+1],posNorm(iInt:iInt+1),'g')
            end
            indx = find(classCurr==2);
            for iInt = indx
                plot([iInt iInt+1],posNorm(iInt:iInt+1),'r')
            end
            indx = find(classCurr==3);
            for iInt = indx
                plot([iInt iInt+1],posNorm(iInt:iInt+1),'c')
            end
            
        end %(if doPlot)
        
    end %(if any(~isnan(protNorm)))
    
end %(for iWin = 1 : numWindows)

%% Edge activity characteristics

%find windows in each category
indxProt  = find(windowMotionType(:,:,1)==1);
indxRet   = find(windowMotionType(:,:,1)==2);
indxPause = find(windowMotionType(:,:,1)==3);
indxUndet = find(isnan(windowMotionType(:,:,1)));

%get fraction of time spent in each category
fracMotionType = [length(indxProt) length(indxRet) length(indxPause) ...
    length(indxUndet)];
fracMotionType = fracMotionType / sum(fracMotionType);

%get median speed in each category
medianSpeedType = [median(avgNormal(indxProt)) median(avgNormal(indxRet)) ...
    median(avgNormal(indxPause)) median(avgNormal(indxUndet))];

%get median motion range for all windows
medianWinMotionRange = nanmedian(motionRange);

%% ~~~ the end ~~~

%% old stuff

%         %assign as significant any intervals that are >= 3 frame steps
%         indx = find(intervalLength >= 3);
%         intervalClass(indx) = sign(intervalProt(indx)); % +1 = protrusion, -1 = retraction
%
%         %go over intervals that are <= 2 frame steps and test for significance
%         indx2 = find(intervalLength<=2);
%         oneMinusPvalue = normcdf(abs(intervalProt(indx2)),0,protVecStd);
%         indx = indx2(oneMinusPvalue>0.95);
%         intervalClass(indx) = sign(intervalProt(indx));% +1 = protrusion, -1 = retraction
%         indx = indx2(oneMinusPvalue<=0.95);
%         intervalClass(indx) = 0; % 0 = pause

%         %determine sign of each protrusion normal
%         protNormSign = sign(protNorm);
%
%         %extract points of switching
%         protNormSignDiff = diff(protNormSign);
%         switchPoints = find(protNormSignDiff~=0) + 1;
%         switchPoints = [1; switchPoints; numFrames]; %#ok<AGROW>
%
%         %calculate length of intervals between switching
%         intervalLength = diff(switchPoints);
%
%         %get number of intervals
%         numIntervals = length(intervalLength);
%
%         %calculate position change (i.e. cumulative protrusion/retraction) in
%         %each interval
%         intervalProt = NaN(numIntervals,1);
%         for iInt = 1 : numIntervals
%             intervalProt(iInt) = sum(protNorm(switchPoints(iInt):switchPoints(iInt+1)-1));
%         end
%
%         %initialize interval classification vector
%         intervalClass = NaN(size(intervalLength));
%
%         %go over all intervals and test for protrusion/retraction significance
%         oneMinusPvalue = normcdf(abs(intervalProt),0,protVecStd);
%         indx = find(oneMinusPvalue>0.9);
%         intervalClass(indx) = sign(intervalProt(indx));% +1 = protrusion, -1 = retraction
%         intervalClass(oneMinusPvalue<=0.9) = 0; % 0 = pause
%
%         %from intervalClass, assign classification to individual frame steps
%         for iInt = 1 : numIntervals
%             classCurr(switchPoints(iInt):switchPoints(iInt+1)-1) = intervalClass(iInt);
%         end
%
%         %         %if the two intervals on either side of a pause interval are both
%         %         %protrusion or both retraction, reclassify the interval to be the
%         %         %same as its immediate neighbors
%         %         indx = find(classCurr(2:end-1)==0)+1;
%         %         indxChange = indx( (classCurr(indx+1) == classCurr(indx-1)) & classCurr(indx+1)~=0 );
%         %         classCurr(indxChange) = classCurr(indxChange-1);

%         %if a tentative pause interval has a protrusion normal > 0.5*protVecStd
%         %and it is surrounded by protrusion on both sides, or if it has a
%         %protrusion normal < -0.5*protVecStd and it is surrounded
%         %by retraction on both sides, reclassify it as protrusion or
%         %retraction
%         indxPause = find(classCurr==0);
%         indxPause = setdiff(indxPause,[1 numFrames-1]);
%         indxProt = indxPause( protNorm(indxPause)'>0.5*protVecStd & ...
%             classCurr(indxPause-1)==1 & classCurr(indxPause+1)==1 );
%         indxRetr = indxPause( protNorm(indxPause)'<-0.5*protVecStd & ...
%             classCurr(indxPause-1)==-1 & classCurr(indxPause+1)==-1 );
%         classCurr(indxProt) = 1;
%         classCurr(indxRetr) = -1;

%         %if the two intervals on either side of a pause interval are both
%         %protrusion or both retraction, reclassify it as a special pause
%         indx = find(classCurr(2:end-1)==0)+1;
%         indxChange = indx( (classCurr(indx+1) == classCurr(indx-1)) & classCurr(indx+1)~=0 );
%         classCurr(indxChange) = classCurr(indxChange-1);
%
%         %if there is a sequence of consecutive pause intervals, but they
%         %all have the same sign, which is the same as the non-pause
%         %intervals to their left and right, then reclassify all the pause
%         %intervals to be the same as their immediate neighbors
%         classCurrTmp = abs(classCurr);
%         classCurrTmp(isnan(classCurrTmp)) = 1;
%         classSwitch = diff(classCurrTmp);
%         pauseStartIndx = find(classSwitch==-1) + 1;
%         pauseEndIndx   = find(classSwitch==1);
%         if classCurr(1) == 0
%             pauseEndIndx = pauseEndIndx(2:end);
%         end
%         if classCurr(end) == 0
%             pauseStartIndx = pauseStartIndx(1:end-1);
%         end
%         pauseDuration = pauseEndIndx - pauseStartIndx + 1;
%         indxKeep = find(pauseDuration>=2);
%         pauseStartIndx = pauseStartIndx(indxKeep);
%         pauseEndIndx = pauseEndIndx(indxKeep);
%         pauseDuration = pauseDuration(indxKeep);
%         numPauses = length(indxKeep);
%         for iPause = 1 : numPauses
%             allSameSign = abs(sum(sign(protNorm(pauseStartIndx(iPause):pauseEndIndx(iPause)))))==pauseDuration(iPause);
%             if allSameSign
%                 pauseSign = sign(protNorm(pauseStartIndx(iPause)));
%                 beforeSign = classCurr(pauseStartIndx(iPause)-1);
%                 afterSign = classCurr(pauseEndIndx(iPause)+1);
%                 if beforeSign==afterSign && pauseSign==beforeSign
%                     classCurr(pauseStartIndx(iPause):pauseEndIndx(iPause)) = beforeSign;
%                 end
%             end
%
%         end

%         for classType = 1 : 3
%             classCurrTmp = (classCurr~=classType);
%             classSwitch = diff(classCurrTmp);
%             classStartIndx = find(classSwitch==-1) + 1;
%             classEndIndx   = find(classSwitch==1);
%             if classCurr(1) == classType
%                 classAft(1:classEndIndx(1)) = classCurr(classEndIndx(1)+1);
%                 durCurr(1:classEndIndx(1))  = classEndIndx(1);
%                 classEndIndx = classEndIndx(2:end);
%             end
%             if classCurr(end) == classType
%                 classBef(classStartIndx(end):end) = classCurr(classStartIndx(end)-1);
%                 durCurr(classStartIndx(end):end)  = numFrames-classStartIndx(end);
%                 classStartIndx = classStartIndx(1:end-1);
%             end
%             numClass = length(classStartIndx);
%             for iClass = 1 : numClass
%                 classBef(classStartIndx(iClass):classEndIndx(iClass)) = classCurr(classStartIndx(iClass)-1);
%                 classAft(classStartIndx(iClass):classEndIndx(iClass)) = classCurr(classEndIndx(iClass)+1);
%                 durCurr(classStartIndx(iClass):classEndIndx(iClass))  = classEndIndx(iClass)-classStartIndx(iClass)+1;
%             end
%         end
        
