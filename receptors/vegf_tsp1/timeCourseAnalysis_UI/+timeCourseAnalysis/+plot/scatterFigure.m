function [ figureHandles ] = scatterFigure( commonInfo, figureData, outputDirFig, doInParallel, outputDirFigA )
%scatterFigure Scatter figure for time course analysis
%
% INPUT
% commonInfo (optional) - structure or string indicating folder of
%                         figureData.mat
% figureData (optional) - structure containing data per figure
%
% outputDirFig (optional) - will pause between figures
%
% If no parameters are given, figureData.mat will be loaded from the
% current directory.
%
% If a single parameter is given, figureData.mat will be loaded from the
% directory indicated.
%
% OUTPUT
% figureHandles - if requested, all figures will stay open
%
% Mark Kittisopikul, November 2015
% based on Tae Kim, Summer 2015

if(nargin < 1)
    [figureData,commonInfo] = timeCourseAnalysis.util.getFigureData;
end

if(nargin < 2 && ischar(commonInfo))
    [figureData,commonInfo] = timeCourseAnalysis.util.getFigureData(commonInfo);
end

if(nargin < 3)
    outputDirFig = [];
end

if(nargin < 4)
    doInParallel = false;
    if(isempty(outputDirFig))
        outputDirFig = pwd;
    end
end

if(nargin < 5)
    outputDirFigA = [];
end

if(nargout > 0)
    noClose = true;
else
    noClose = false;
end

colors = timeCourseAnalysis.plot.getColors(commonInfo.conditions);

if(~isfield(figureData,'inOutFlag'))
    [figureData.inOutFlag] = deal([]);
end

if ~exist(outputDirFig, 'dir')
    mkdir(outputDirFig);
end
if ~exist(outputDirFigA, 'dir')
    mkdir(outputDirFigA);
end

% doInParallel = 0;
if(doInParallel)
    figureHandles = pararrayfun_progress( ...
        @(f) scatterIndividualFigure( ...
        commonInfo, ...
        f, ...
        colors, ...
        outputDirFig, ...
        noClose, ...
        outputDirFigA), ...
        figureData, ...
        'UniformOutput', false, ...
        'ErrorHandler',@scatterFigureErrorParallel ...
        ,'UseErrorStruct',false ...
        );
else
    figureHandles = arrayfun( ...
        @(f) scatterIndividualFigure( ...
        commonInfo, ...
        f, ...
        colors, ...
        outputDirFig, ...
        noClose, ...
        outputDirFigA), ...
        figureData, ...
        'UniformOutput', false, ...
        'ErrorHandler',@scatterFigureError ...
        );
end

if(nargout > 0)
    figureHandles = [figureHandles{:}];
end

end
function figureHandle = scatterIndividualFigure(commonInfo, figureData, colors, outputDirFig, noClose, outputDirFigA)

plotTitle = [figureData.titleBase ' ' figureData.titleVariable];


%% KJ: This is Tae's original spline fit (with some modifications)

figureHandle = figure('Name', plotTitle);
hold on;

%plot all inlier data and store line handles
inIdx = cellfun(@(inOutFlag) inOutFlag==1,figureData.inOutFlag,'UniformOutput',false);
lineHandleP1 = cellfun(@(fitData,times,data,inIdx) plot(fitData,times(inIdx),data(inIdx)),figureData.fitData,commonInfo.times,figureData.data,inIdx,'UniformOutput',false);

%KJ: plot outliers not used in fit
if(~isempty(figureData.inOutFlag))
    %outliers based on value
    outIdx = cellfun(@(inOutFlag) inOutFlag==0,figureData.inOutFlag,'UniformOutput',false);
    lineHandle0 = cellfun(@(times,data,outIdx) plot(times(outIdx),data(outIdx),'o'),commonInfo.times,figureData.data,outIdx,'UniformOutput',false);
    %outliers based on isolation
    outIdx = cellfun(@(inOutFlag) inOutFlag==-1,figureData.inOutFlag,'UniformOutput',false);
    lineHandleM1 = cellfun(@(times,data,outIdx) plot(times(outIdx),data(outIdx),'s'),commonInfo.times,figureData.data,outIdx,'UniformOutput',false);
end

%plot vertical lines indicating aligning Times
nCond = numel(figureData.data);
for iCond = 1:nCond
    plot([commonInfo.timeShift(iCond), commonInfo.timeShift(iCond)], [figureData.yMax, figureData.yMin], 'Color', colors{iCond});
end

%change the color so that color of data and fit match
cellfun(@(x, y) set(x, 'Color', y), lineHandleP1, colors);
cellfun(@(x, y) set(x, 'Color', y), lineHandle0, colors);
cellfun(@(x, y) set(x, 'Color', y), lineHandleM1, colors);

%create legends
fitHandle = [lineHandleP1{:}];
legend(fitHandle(2,:), commonInfo.conditions);

%axis limit
if ~isempty(figureData.yMax)
    ylim([figureData.yMin, figureData.yMax]);
end

%label axis
xlabel('Time (min)');
ylabel(figureData.yLabel);
title(plotTitle);

%save and close
if(~isempty(outputDirFig))
    savefig(figureHandle, [outputDirFig filesep plotTitle '.fig']);
else
    if(~noClose)
        pause;
    end
end
if(~noClose)
    close(figureHandle);
end

%% KJ: New "moving average" analysis

figureHandle = figure('Name', plotTitle);
hold on;

%collect all moving average information
for i=1:length(figureData.dataAve)
    timesTmp{i,1} = figureData.dataAve{i}.timeAve;
    dataTmp{i,1} = figureData.dataAve{i}.dataAve;
    semTmp = figureData.dataAve{i}.dataStd ./ sqrt(figureData.dataAve{i}.nDataPts);
    dataPlusTmp{i,1} = dataTmp{i} + semTmp;
    dataMinusTmp{i,1} = dataTmp{i} - semTmp;
end

%plot all inlier data and store line handles
inIdx = cellfun(@(inOutFlag) inOutFlag==1,figureData.inOutFlag,'UniformOutput',false);
lineHandleP1 = cellfun(@(times,data,inIdx) plot(times(inIdx),data(inIdx),'.'), commonInfo.times, figureData.data, inIdx, 'UniformOutput', false);
lineHandleA  = cellfun(@(times,data) plot(times,data),timesTmp,dataTmp,'UniformOutput',false);
lineHandleAP = cellfun(@(times,data) plot(times,data,'LineStyle',':'),timesTmp,dataPlusTmp,'UniformOutput',false);
lineHandleAM = cellfun(@(times,data) plot(times,data,'LineStyle',':'),timesTmp,dataMinusTmp,'UniformOutput',false);

%KJ: plot outliers not used in fit
if(~isempty(figureData.inOutFlag))
    %outliers based on value
    outIdx = cellfun(@(inOutFlag) inOutFlag==0,figureData.inOutFlag,'UniformOutput',false);
    lineHandle0 = cellfun(@(times,data,outIdx) plot(times(outIdx),data(outIdx),'o'),commonInfo.times,figureData.data,outIdx,'UniformOutput',false);
    %outliers based on isolation
    outIdx = cellfun(@(inOutFlag) inOutFlag==-1,figureData.inOutFlag,'UniformOutput',false);
    lineHandleM1 = cellfun(@(times,data,outIdx) plot(times(outIdx),data(outIdx),'s'),commonInfo.times,figureData.data,outIdx,'UniformOutput',false);
end

%plot vertical lines indicating aligning Times
nCond = numel(figureData.data);
for iCond = 1:nCond
    plot([commonInfo.timeShift(iCond), commonInfo.timeShift(iCond)], [figureData.yMax, figureData.yMin], 'Color', colors{iCond});
end

%change the color so that color of data and fit match
cellfun(@(x, y) set(x, 'Color', y), lineHandleP1, colors);
cellfun(@(x, y) set(x, 'Color', y), lineHandle0, colors);
cellfun(@(x, y) set(x, 'Color', y), lineHandleM1, colors);
cellfun(@(x, y) set(x, 'Color', y), lineHandleA, colors);
cellfun(@(x, y) set(x, 'Color', y), lineHandleAP, colors);
cellfun(@(x, y) set(x, 'Color', y), lineHandleAM, colors);

%create legends
fitHandle = [lineHandleA{:}];
legend(fitHandle, commonInfo.conditions);

%axis limit
if ~isempty(figureData.yMax)
    ylim([figureData.yMin, figureData.yMax]);
end

%label axis
xlabel('Time (min)');
ylabel(figureData.yLabel);
title(plotTitle);

%save and close
if(~isempty(outputDirFigA))
    savefig(figureHandle, [outputDirFigA filesep plotTitle '.fig']);
else
    if(~noClose)
        pause;
    end
end
if(~noClose)
    close(figureHandle);
end

end

function out = scatterFigureError(err,figureData)
if(isempty(err.identifier))
    error('scatterFigure:scatterFigureError',['Could not plot ' figureData.titleBase ' ' figureData.titleVariable])
end
try
    warning('scatterFigure:scatterFigureError',['Could not plot ' figureData.titleBase ' ' figureData.titleVariable]);
    disp(err);
    disp(figureData);
    out = [];
    disp('Continuing ...');
catch newError
    disp('scatterFigureError:ErrorInErrorFunction','An error has occured in timeCourseAnalysis.scatterFigure/scatterFigureError');
    disp(getReport(newError));
end
end

function out = scatterFigureErrorParallel(err,d)
% Take advantage of the newer MException class to print more
% informative errors
try
    figureData = d.InputArguments{1};
    warning('scatterFigure:scatterFigureError',['Could not plot ' figureData.titleBase ' ' figureData.titleVariable]);
    disp(getReport(err));
    disp(d);
    disp(figureData);
    out = [];
    disp('Continuing ...');
catch newError
    warning('scatterFigureErrorParallel:ErrorInErrorFunction','An error has occured in timeCourseAnalysis.scatterFigure/scatterFigureErrorParallel');
    try
        disp(getReport(newError));
    catch newNewError
        warning('scatterFigureErrorParallel:ErrorInErrorFunction2','A deeper error has occured in timeCourseAnalysis.scatterFigure/scatterFigureErrorParallel');
    end
end
end
