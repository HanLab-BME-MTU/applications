function [fracMotionType,motionTypeChar,winMotionRange,autoCorrProtComb,...
    autoCorrProtInd,indxGood] = characterizeEdgeMotion(avgNormal,...
    windowMotionType,doPlot,savePlotDir)
%CHARACTERIZEEDGEMOTION characterizes edge motion in each window based on protrusion vector samples
%
%SYNPOSIS [fracMotionType,motionTypeChar,winMotionRange,autoCorrProtComb,...
%    autoCorrProtInd] = characterizeEdgeMotion(avgNormal,windowMotionType,doPlot)
%
%INPUT  avgNormal       : The protrusion vector normals for each window.
%       windowMotionType: Output of classifyEdgeMotion.
%       doPlot          : Flag with value 1 to plot speed and duration
%                         histograms and 0 to not plot.
%                         Optional. Default: 0.
%       savePlotDir     : Directory where plots are to be saved.
%                         Optional. If not supplied, plots will not be
%                         saved.
%OUTPUT fracMotionType  : Fraction of intervals classified as protrusion,
%                         retraction, pause and undetermined. 
%                         1st row: all windows.
%                         2nd row: only "good" windows.
%       motionTypeChar  : Structure with 2 fields, "all" and "good",
%                         with results from all windows and only good
%                         windows. "all" and "good" are structure arrays
%                         with 3 entries for protrusion, retraction and
%                         pause, each with the fields:
%           .speed          : Speed distribution.
%           .duration       : Duration distribution (i.e. time from onset
%                             until switching to another motion type).
%       winMotionRange  : Motion range of each window.
%       autoCorrProtComb: Structure with 2 fields, "all" and "good",
%                         with results from all windows and only good
%                         windows. "all"/"good" are arrays storing the
%                         protrusion normal autocorrelation for all/good
%                         windows combined.
%       autoCorrProtInd : Protrusion normal autocorrelation for individual
%                         windows (each window is a layer in the 3rd
%                         dimension).
%       indxGood        : Indices of good windows, defined as windows whose
%                         protrusion normal autocorrelation is NOT negative
%                         at lag 1 and then positive at lag 2.
%
%Khuloud Jaqaman, July 2012

%% Input

if nargin < 2
    error('--characterizeEdgeMotion: Please enter protrusion samples and window activity classification.');
end

if nargin < 3 || isempty(doPlot)
    doPlot = 0;
end

if nargin < 4 || isempty(savePlotDir)
    savePlotDir = [];
end

%get number of windows and number of frames
[numWindows,numFrames] = size(avgNormal);

%% Identify good windows based on individual window autocorrelation

%get maximum time lag for autocorrelation
maxLagInd = floor(numFrames/4);
maxLagCombMax = max(floor(numFrames/2),10);
maxLagCombMin = min(floor(numFrames/2),10);

%calculate protrusion normal autocorrelation for individual windows
protNormAllWindows = repmat(struct('observations',[]),numWindows,1);
autoCorrProtInd = NaN(maxLagInd+1,2,numWindows);
for iWindow = 1 : numWindows
    protNormAllWindows(iWindow).observations = avgNormal(iWindow,:)';
    autoCorrProtInd(:,:,iWindow) = autoCorr(protNormAllWindows(iWindow),maxLagInd);
end

%identify windows with "good" protrusion activity, defined as those
%windows whose autocorrelation is NOT negative at lag 1 and then positive at lag 2
if size(autoCorrProtInd,1) >= 3
    autoCorrLags12 = squeeze(autoCorrProtInd(2:3,1,:))';
    indxBad = find( (autoCorrLags12(:,1)<0 & autoCorrLags12(:,2)>0) | isnan(autoCorrLags12(:,1)) );
else
    indxBad = (1:numWindows)';
end
indxGood = setdiff((1:numWindows)',indxBad);

%% Combined autocorrelation

%calculate protrusion normal autocorrelation for all windows combined
[autoCorrProtCombAll,maxLagCombUsedAll] = autoCorr3attempts(protNormAllWindows,maxLagCombMax,maxLagCombMin,maxLagInd);

%calculate the combined autocorrelation but only for "good" windows
if ~isempty(indxGood)
    [autoCorrProtCombGood,maxLagCombUsedGood] = autoCorr3attempts(protNormAllWindows(indxGood),maxLagCombMax,maxLagCombMin,maxLagInd);
else
    autoCorrProtCombGood = NaN(1,2);
    maxLagCombUsedGood = 1;
end

autoCorrProtComb.all = autoCorrProtCombAll;
autoCorrProtComb.good = autoCorrProtCombGood;

%% Characteristics from activity classification

%get activity characteristics for all windows
[motionTypeCharAll,fracMotionType(1,:)] = getCharLocal(windowMotionType,avgNormal);

%get activity characteristics for good windows
if ~isempty(indxGood)
    [motionTypeCharGood,fracMotionType(2,:)] = getCharLocal(windowMotionType(indxGood,:,:),avgNormal(indxGood,:));
else
    motionTypeCharGood = [];
    fracMotionType(2,:) = NaN;
end

%put output together
motionTypeChar.all = motionTypeCharAll;
motionTypeChar.good = motionTypeCharGood;

%get motion range for each window
winMotionRange = NaN(numWindows,1);
for iWindow = 1 : numWindows
    protNorm = avgNormal(iWindow,1:end-1)';
    if any(~isnan(protNorm))
        posNorm = nancumsumwithzero(protNorm);
        winMotionRange(iWindow) = max(posNorm)-min(posNorm);
    end
end

%% Plotting

if doPlot

    hFig = figure('Name','Speeds');
    subplot(1,2,1)
    hold on
    speedProt = vertcat(motionTypeCharAll(1).speed);
    speedRet = vertcat(motionTypeCharAll(2).speed);
    speedPause = vertcat(motionTypeCharAll(3).speed);
    dataLength = [length(speedProt) length(speedRet) length(speedPause)];
    dataAll = NaN(max(dataLength),3);
    dataAll(1:dataLength(1),1) = speedProt;
    dataAll(1:dataLength(2),2) = speedRet;
    dataAll(1:dataLength(3),3) = speedPause;
    divStep = (max(dataAll(:))-min(dataAll(:)))/30;
    binEdge = min(dataAll(:)):divStep:max(dataAll(:));
    binCenter = mean([binEdge(1:end-1);binEdge(2:end)]);
    distributionPlot(dataAll,'color',[0.9 0.9 0.9],'histOpt',0,'divFactor',binCenter,'showMM',2)
    boxplot(dataAll,'notch','on','symbol','g.')
    title('All windows')
    hold off
    if ~isempty(indxGood)
        subplot(1,2,2)
        hold on
        speedProt = vertcat(motionTypeCharGood(1).speed);
        speedRet = vertcat(motionTypeCharGood(2).speed);
        speedPause = vertcat(motionTypeCharGood(3).speed);
        dataLength = [length(speedProt) length(speedRet) length(speedPause)];
        dataAll = NaN(max(dataLength),3);
        dataAll(1:dataLength(1),1) = speedProt;
        dataAll(1:dataLength(2),2) = speedRet;
        dataAll(1:dataLength(3),3) = speedPause;
        divStep = (max(dataAll(:))-min(dataAll(:)))/30;
        binEdge = min(dataAll(:)):divStep:max(dataAll(:));
        binCenter = mean([binEdge(1:end-1);binEdge(2:end)]);
        distributionPlot(dataAll,'color',[0.9 0.9 0.9],'histOpt',0,'divFactor',binCenter,'showMM',2)
        boxplot(dataAll,'notch','on','symbol','g.')
        title('Good windows')
        hold off
    end
    if ~isempty(savePlotDir)
        saveas(hFig,fullfile(savePlotDir,'activitySpeed.fig'));
    end
    
    hFig = figure('Name','Durations');
    subplot(1,2,1)
    hold on
    durationProt = vertcat(motionTypeCharAll(1).duration);
    durationRet = vertcat(motionTypeCharAll(2).duration);
    durationPause = vertcat(motionTypeCharAll(3).duration);
    dataLength = [length(durationProt) length(durationRet) length(durationPause)];
    dataAll = NaN(max(dataLength),3);
    dataAll(1:dataLength(1),1) = durationProt;
    dataAll(1:dataLength(2),2) = durationRet;
    dataAll(1:dataLength(3),3) = durationPause;
    distributionPlot(dataAll,'color',[0.9 0.9 0.9],'histOpt',0,'divFactor',1:numFrames,'showMM',2)
    boxplot(dataAll,'notch','on','symbol','g.')
    title('All windows')
    hold off
    if ~isempty(indxGood)
        subplot(1,2,2)
        hold on
        durationProt = vertcat(motionTypeCharGood(1).duration);
        durationRet = vertcat(motionTypeCharGood(2).duration);
        durationPause = vertcat(motionTypeCharGood(3).duration);
        dataLength = [length(durationProt) length(durationRet) length(durationPause)];
        dataAll = NaN(max(dataLength),3);
        dataAll(1:dataLength(1),1) = durationProt;
        dataAll(1:dataLength(2),2) = durationRet;
        dataAll(1:dataLength(3),3) = durationPause;
        distributionPlot(dataAll,'color',[0.9 0.9 0.9],'histOpt',0,'divFactor',1:numFrames,'showMM',2)
        boxplot(dataAll,'notch','on','symbol','g.')
        title('Good windows')
        hold off
    end
    if ~isempty(savePlotDir)
        saveas(hFig,fullfile(savePlotDir,'activityDuration.fig'));
    end
    
    hFig = figure('Name','Protrusion normal autocorrelation - combined');
    subplot(1,2,1)
    hold on
    plot(0:maxLagCombUsedAll,autoCorrProtCombAll(:,1))
    myErrorbar(0:maxLagCombUsedAll,autoCorrProtCombAll(:,1),autoCorrProtCombAll(:,2))
    tmp = vertcat(protNormAllWindows.observations);
    numPointsAuto = length(find(~isnan(tmp)));
    plot([-0.1 maxLagCombUsedAll+0.1],1.96/sqrt(numPointsAuto)*[1 1],'k')
    plot([-0.1 maxLagCombUsedAll+0.1],-1.96/sqrt(numPointsAuto)*[1 1],'k')
    xlim([-0.1 maxLagCombUsedAll+0.1])
    ylim([-1.1 1.1])
    xlabel('Time lag (frames)')
    ylabel('Auto-correlation')
    title('All windows')
    hold off
    if ~isempty(indxGood)
        subplot(1,2,2)
        hold on
        plot(0:maxLagCombUsedGood,autoCorrProtCombGood(:,1))
        myErrorbar(0:maxLagCombUsedGood,autoCorrProtCombGood(:,1),autoCorrProtCombGood(:,2))
        tmp = vertcat(protNormAllWindows(indxGood).observations);
        numPointsAuto = length(find(~isnan(tmp)));
        plot([-0.1 maxLagCombUsedGood+0.1],1.96/sqrt(numPointsAuto)*[1 1],'k')
        plot([-0.1 maxLagCombUsedGood+0.1],-1.96/sqrt(numPointsAuto)*[1 1],'k')
        xlim([-0.1 maxLagCombUsedGood+0.1])
        ylim([-1.1 1.1])
        xlabel('Time lag (frames)')
        ylabel('Auto-correlation')
        title('Good windows')
        hold off
    end
    if ~isempty(savePlotDir)
        saveas(hFig,fullfile(savePlotDir,'protAutocorrComb.fig'));
    end
    
    hFig = figure('Name','Protrusion normal autocorrelation - individual');
    subplot(1,2,1)
    hold on
    cmapExt = [[0 0 0]; colormap];
    autoCorr2plot = squeeze(autoCorrProtInd(:,1,:))';
    nanValue = -1 - 2/64;
    autoCorr2plot(isnan(autoCorr2plot)) = nanValue;
    imagesc((0:maxLagInd),(1:numWindows),autoCorr2plot);
    caxis([nanValue 1])
    colormap(cmapExt);
    colorbar
    xlabel('Time lag (frames)')
    ylabel('Position along edge (window #)')
    hold off
    subplot(1,2,2)
    hold on
    if ~isempty(indxGood)
        autoCorr2plotGood = squeeze(autoCorrProtInd(:,1,indxGood));
        plot(0:maxLagInd,autoCorr2plotGood(:,1),'r')
    end
    autoCorr2plotBad = squeeze(autoCorrProtInd(:,1,indxBad));
    plot(0:maxLagInd,autoCorr2plotBad(:,1),'b')
    if ~isempty(indxGood)
        legend({'Good','Bad'})
    else
        legend({'Bad'})
    end
    plot(0:maxLagInd,autoCorr2plotBad,'b')
    if ~isempty(indxGood)
        plot(0:maxLagInd,autoCorr2plotGood,'r')
    end
    xlim([-0.1 maxLagInd+0.1])
    ylim([-1.1 1.1])
    xlabel('Time lag (frames)')
    ylabel('Auto-correlation')    
    hold off
    if ~isempty(savePlotDir)
        saveas(hFig,fullfile(savePlotDir,'protAutocorrInd.fig'));
    end
    
end

%% Sub-function 1

function [autoCorrProtComb,maxLagCombUsed] = autoCorr3attempts(protNormAllWindows,maxLagCombMax,maxLagCombMin,maxLagInd)

try
    autoCorrProtComb = autoCorr(protNormAllWindows,maxLagCombMax,-1,1);
catch
    try
        autoCorrProtComb = autoCorr(protNormAllWindows,maxLagCombMin,-1,1);
    catch
        autoCorrProtComb = autoCorr(protNormAllWindows,maxLagInd,-1,1);
    end
end
maxLagCombUsed = size(autoCorrProtComb,1) - 1;


%% Sub-function 2

function [motionTypeChar,fracMotionType] = getCharLocal(windowMotionType,avgNormal)

%get number of windows and number of frames
[numWindows,numFrames] = size(avgNormal);

%find windows in each activity category
indxProt  = find(windowMotionType(:,:,1)==1);
indxRet   = find(windowMotionType(:,:,1)==2);
indxPause = find(windowMotionType(:,:,1)==3);
indxUndet = find(isnan(windowMotionType(:,:,1)));

%get fraction of time spent in each activity
fracMotionType = [length(indxProt) length(indxRet) length(indxPause) ...
    length(indxUndet)];
fracMotionType = fracMotionType / sum(fracMotionType);

%get speed distribution in each activity
speedProt = avgNormal(indxProt);
speedRet = avgNormal(indxRet);
speedPause = avgNormal(indxPause);

%get duration distribution of each activity
%ignore events that are there from the beginning of the movie or until the
%end of the movie, as it is unknown when they started or when they ended
classDur = repmat(struct('values',[]),3,1);
for iWindow = 1 : numWindows
    classCurr = windowMotionType(iWindow,:,1);
    for iClass = 1 : 3
        classTrans = diff([0 classCurr==iClass 0]);
        indxClassS = find(classTrans==1);
        indxClassE = find(classTrans==-1)-1;
        indxKeep = find( (indxClassS~=1) & (indxClassE~=(numFrames-1)) );
        indxClassS = indxClassS(indxKeep); %#ok<FNDSB>
        %         indxClassE = indxClassE(indxKeep);
        classDur(iClass).values = [classDur(iClass).values; windowMotionType(iWindow,indxClassS,4)'];
    end
end
durationProt = vertcat(classDur(1).values);
durationRet = vertcat(classDur(2).values);
durationPause = vertcat(classDur(3).values);

%assemble characteristics for output
motionTypeChar(1).speed = speedProt;
motionTypeChar(1).duration = durationProt;
motionTypeChar(2).speed = speedRet;
motionTypeChar(2).duration = durationRet;
motionTypeChar(3).speed = speedPause;
motionTypeChar(3).duration = durationPause;


%% ~~~ the end ~~~

