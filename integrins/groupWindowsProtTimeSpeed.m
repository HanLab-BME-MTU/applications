function windowGroup = groupWindowsProtTimeSpeed(protSamples,doPlot,...
    windowRange,frameRange)
%CLASSIFYEDGEMOTION classifies edge motion in each window based on protrusion vector samples
%
%SYNPOSIS windowGroup = groupWindowsProtTimeSpeed(protSamples,doPlot,...
%    windowRange,frameRange)
%
%INPUT  protSample      : The protrusion samples as output by the windowing
%                         software.
%       doPlot          : Flag with value 1 to plot and 0 not to plot
%                         classification. Optional. Default: 0.
%       windowRange     : 2-row vector indicating range of windows (i.e.
%                         window number along cell edge) to include in
%                         analysis.
%                         Optional. Default: [].
%       frameRange      : 2-row vector indicating range of window frames
%                         to include in analysis.
%                         Optional. Default: [].
%OUTPUT windowGroup     : 
%
%Khuloud Jaqaman, January 2012

%% Input

if nargin < 1
    error('--groupWindowsProtTimeSpeed: Please enter protrusion samples!');
end

if nargin < 2 || isempty(doPlot)
    doPlot = 0;
end

if nargin < 3 || isempty(windowRange)
    windowRange = [];
end

if nargin < 4 || isempty(frameRange)
    frameRange = [];
end

%% Grouping

%classify the activity of each window
windowMotionType = classifyEdgeMotion(protSamples,doPlot,windowRange,frameRange);

%get number of windows and frames
[numWindows,numFrames] = size(windowMotionType);

%add a column of NaN to windowMotionType to make subsequent steps simpler
windowMotionType = [NaN(numWindows,1) windowMotionType];

%initialize output
windowGroup = repmat(struct('frameWinIndx',NaN(0,2)),numFrames-1,numFrames-1);

%group windows based on protrusion duration
for iDuration = 1 : numFrames-1
    
    %sum up edge motion classification values
    motionTypeSum = windowMotionType(:,2:end-iDuration);
    for iInc = 1 : iDuration-1
        motionTypeSum = motionTypeSum + windowMotionType(:,2+iInc:end-iDuration+iInc);
    end
    
    %find protrusions that last for iDuration intervals
    [winIndxTmp,frameIndxTmp] = find( motionTypeSum==iDuration & ...
        windowMotionType(:,1:end-iDuration-1)~=1 & ...
        windowMotionType(:,iDuration+2:end)~=1 );
    
    %go back in time and sub-classify windows based on events before
    %protrusion onset
    for iTimeBack = 1 : numFrames-1
        
        %find all windows with protrusion onset that does not allow going
        %back enough in time
        indxRemove = find(frameIndxTmp==iTimeBack);
        indxKeep = setdiff((1:length(frameIndxTmp))',indxRemove);
        winIndxTmpGroup1 = winIndxTmp(indxRemove);
        frameIndxTmpGroup1 = frameIndxTmp(indxRemove);
        winIndxTmp = winIndxTmp(indxKeep);
        frameIndxTmp = frameIndxTmp(indxKeep);
        
        %for the remaining windows, look iTimeBack frames and check whether
        %the event is protrusion
        linearIndx = sub2ind([numWindows numFrames],winIndxTmp,frameIndxTmp-iTimeBack+1);
        windowMotionTypeInPast = windowMotionType(linearIndx);
        indxRemove = find(windowMotionTypeInPast==1);
        indxKeep = setdiff((1:length(frameIndxTmp))',indxRemove);
        winIndxTmpGroup1 = [winIndxTmpGroup1; winIndxTmp(indxRemove)]; %#ok<AGROW>
        frameIndxTmpGroup1 = [frameIndxTmpGroup1; frameIndxTmp(indxRemove)]; %#ok<AGROW>
        winIndxTmp = winIndxTmp(indxKeep);
        frameIndxTmp = frameIndxTmp(indxKeep);

        if ~isempty(frameIndxTmpGroup1)
            windowGroup(iDuration,iTimeBack).frameWinIndx = [frameIndxTmpGroup1 winIndxTmpGroup1];
        end
        
    end
    
end



%% ~~~ the end ~~~

 