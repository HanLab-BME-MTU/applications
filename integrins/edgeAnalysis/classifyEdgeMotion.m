function [windowMotionType,windowMotionChar] = classifyEdgeMotion(...
    protSamples,doPlot,indxWindows,frameRange,scheme,savePlotDir)
%CLASSIFYEDGEMOTION classifies edge motion in each window based on protrusion vector samples
%
%SYNPOSIS [windowMotionType,windowMotionChar] = classifyEdgeMotion(...
%    protSamples,doPlot,indxWindows,frameRange,scheme,savePlotDir)
%
%INPUT  protSample      : The protrusion samples as output by the windowing
%                         software.
%       doPlot          : Flag with values from -1 to 2:
%                         0 : plot nothing.
%                         1 : plot invidual windows classification.
%                         2 : 1 + plot window activity characteristics.
%                         -1: plot only window activity characteristics. 
%                         Optional. Default: 0.
%       indxWindows     : Vector with indices of windows (i.e. window number 
%                         along cell edge) to include in analysis.
%                         Optional. Default: [1 (last window)].
%       frameRange      : 2-row vector indicating range of window frames
%                         to include in analysis.
%                         Optional. Default: [1 (last frame)].
%       scheme          : Classification scheme:
%                         ** 0 : most detailed, with all pauses retained.
%                         ** 1 : a bit less detailed, where pauses
%                         surrounded by periods of the same type are
%                         converted to that type.
%                         ** 2 : 1 + protrusion focused, i.e. intervals are
%                         classified as protrusion or not protrusion,
%                         lumping retraction and pause together.
%                         ** 3 : 1 + retraction focused, i.e. intervals are
%                         classified as retraction or not retractions,
%                         lumping protrusion and pause together.
%                         Optional. Default: 0.
%       savePlotDir     : Directory where plots are to be saved. 
%                         Note: Only plots of activity characteristics are
%                         saved.
%                         Optional. If not supplied, plots will not be
%                         saved.
%OUTPUT windowMotionType: Size(protSamples.avgNormal)-by-6 array
%                         indicating motion classification, where 1 means
%                         protrusion, 2 means retraction, 3 means pause,
%                         and NaN indicates windows without a prutrusion
%                         vector.
%                         Classification scheme 2: both retraction and
%                         pause are lumped together as "2".
%                         Classification scheme 3: both protrusion and
%                         pause are lumped together as "1".
%                         In the 3rd dimension:
%                         Entry 1: classification of each interval.
%                         Entry 2: classification of previous period.
%                         Note: period = group of consecutive intervals
%                         with same classification.
%                         Entry 3: classification of next period.
%                         Entry 4: duration of the period that the interval
%                         belongs to.
%                         Entry 5: duration of previous period.
%                         Entry 6: duration of next period.
%       windowMotionChar: Structure with fields:
%           .fracMotionType  : Fraction of intervals classified as protrusion,
%                              retraction, pause and undetermined. 
%                              1st row: all windows.
%                              2nd row: only "good" windows.
%           .motionTypeChar  : Structure with 2 fields, "all" and "good",
%                              with results from all windows and only good
%                              windows. "all" and "good" are structure arrays
%                              with 3 entries for protrusion, retraction and
%                              pause, each with the fields:
%               .speed          : Speed distribution.
%               .duration       : Duration distribution (i.e. time from onset
%                                 until switching to another motion type).
%           .winMotionRange  : Motion range of each window.
%           .autoCorrProtComb: Structure with 2 fields, "all" and "good",
%                              with results from all windows and only good
%                              windows. "all"/"good" are arrays storing the
%                              protrusion normal autocorrelation for all/good
%                              windows combined.
%           .autoCorrProtInd : Protrusion normal autocorrelation for individual
%                              windows (each window is a layer in the 3rd
%                              dimension).
%           .indxGood        : Indices of good windows, defined as windows whose
%                              protrusion normal autocorrelation is NOT negative
%                              at lag 1 and then positive at lag 2.
%
%
%Khuloud Jaqaman, December 2011

%% Input

if nargin < 1
    error('--classifyEdgeMotion: Please enter protrusion samples!');
end

if nargin < 2 || isempty(doPlot)
    doPlot = 0;
end

%get number of windows and number of frames
avgNormal = protSamples.avgNormal;
[numWindows,numFrames] = size(avgNormal);

if nargin < 3 || isempty(indxWindows)
    indxWindows = (1:numWindows)';
end

if nargin < 4 || isempty(frameRange)
    frameRange = [1 numFrames];
end

if nargin < 5 || isempty(scheme)
    scheme = 0;
end

if nargin < 6 || isempty(savePlotDir)
    savePlotDir = [];
end

%% Classification

%for windows/frames outside of the range of interest, convert their values
%to NaN
indxWinRemove = setdiff((1:numWindows)',indxWindows);
avgNormal(indxWinRemove,:,:) = NaN;
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
        
        %modify classification based on classification scheme
        if scheme > 0
            
            %convert pauses of length 1 frame and surrounded by protrusion
            %on both sides or retraction on both sides to the
            %classification of their neighbors
            pauseTrans = diff([0 classCurr==3 0]);
            indxPauseS = find(pauseTrans==1);
            indxPauseE = find(pauseTrans==-1)-1;
            %             indxKeep = find( (indxPauseS~=1) & (indxPauseE~=(numFrames-1)) );
            indxKeep = find( (indxPauseS~=1) & (indxPauseE~=(numFrames-1)) & ((indxPauseE-indxPauseS)==0));
            indxPauseS = indxPauseS(indxKeep);
            indxPauseE = indxPauseE(indxKeep);
            indxConvert = find(classCurr(indxPauseS-1)==classCurr(indxPauseE+1));
            for iPause = indxConvert
                classCurr(indxPauseS(iPause):indxPauseE(iPause)) = classCurr(indxPauseS(iPause)-1);
            end
            
            switch scheme
                case 2
                    classCurr(classCurr==3) = 2;
                case 3
                    classCurr(classCurr==3) = 1;
            end
            
        end
        
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
        
        if (doPlot == 1 || doPlot == 2)
            
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
            hold off
            
        end %(if doPlot)
        
    end %(if any(~isnan(protNorm)))
    
end %(for iWin = 1 : numWindows)

%% Edge activity characteristics

if (doPlot == -1 || doPlot == 2)
    doPlot2 = 1;
else
    doPlot2 = 0;
end
[fracMotionType,motionTypeChar,winMotionRange,autoCorrProtComb,...
    autoCorrProtInd,indxGoodWindow] = characterizeEdgeMotion(avgNormal,...
    windowMotionType,doPlot2,savePlotDir);

windowMotionChar = struct('fracMotionType',fracMotionType,...
    'motionTypeChar',motionTypeChar,'winMotionRange',winMotionRange,...
    'autoCorrProtComb',autoCorrProtComb,'autoCorrProtInt',autoCorrProtInd,...
    'indxGoodWindow',indxGoodWindow);


%% ~~~ the end ~~~

