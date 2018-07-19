function [ figureHandles ] = standardErrorFigure( commonInfo, figureData, scatterData, outputDirFig2, doInParallel )
%standardErrorFigure Plots time course analysis figure with standard error
%
% INPUT
% commonInfo (optional) - structure or string indicating folder of
%                         figureData.mat
% figureData (optional) - structure containing data per figure
%
% outputDirFig (optional) - will pause between figures
% scatterData - logical whether to scatter data (default: true)
% outputDirFig2 - directory to save the figure
%
% If no parameters are given, figureData.mat will be loaded from the
% current directory.
%
% If a single parameter is given, figureData.mat will be loaded from the
% directory indicated.
%
% OUTPUT
% none
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
    scatterData = true;
end
if(nargin < 4)
    outputDirFig2 = [];
end
if(nargin < 5)
    doInParallel = false;
    if(isempty(outputDirFig2))
        outputDirFig2 = pwd;
    end
end

if ~exist(outputDirFig2, 'dir')
    mkdir(outputDirFig2);
end

if(nargout > 0)
    noClose = true;
else
    noClose = false;
end

colors = timeCourseAnalysis.plot.getColors(commonInfo.conditions);

% doInParallel = 0;
if(doInParallel)
    figureHandles = pararrayfun_progress( ...
        @(f) standardErrorIndividualFigure( ...
        commonInfo, ...
        f, ...
        colors, ...
        scatterData, ...
        outputDirFig2, ...
        noClose), ...
        figureData, ...
        'UniformOutput', false, ...
        'ErrorHandler',@standardErrorFigureErrorParallel ...
        ,'UseErrorStruct',false ...
        );
else
    figureHandles = arrayfun( ...
        @(f) standardErrorIndividualFigure( ...
        commonInfo, ...
        f, ...
        colors, ...
        scatterData, ...
        outputDirFig2, ...
        noClose), ...
        figureData, ...
        'UniformOutput', false, ...
        'ErrorHandler',@standardErrorFigureError ...
        );
end

if(nargout > 0)
    figureHandles = [figureHandles{:}];
end

end
function figureHandle = standardErrorIndividualFigure(commonInfo, figureData, colors, scatterData, outputDirFig2,noClose)

%open up the figure and make it current figure
nConditions = numel(commonInfo.conditions);
lineObj = cell(1,nConditions);
figureHandle = figure();
hold on;

%plot standard errors
for iCond = 1:nConditions
    
    fitValues = figureData.fitData{iCond}(commonInfo.analysisTimes{iCond});
    fitValues = fitValues';
    
    %KJ: plot inlier data
    plot(commonInfo.analysisTimes{iCond}, fitValues + figureData.fitError{iCond}, ':', 'color', colors{iCond});
    plot(commonInfo.analysisTimes{iCond}, fitValues - figureData.fitError{iCond}, ':', 'color', colors{iCond});
    lineObj{iCond} = plot(commonInfo.analysisTimes{iCond}, fitValues, 'color', colors{iCond});
    plot([commonInfo.timeShift(iCond), commonInfo.timeShift(iCond)], [figureData.yMax, figureData.yMin], 'Color', colors{iCond});
    
    inIdx = figureData.inOutFlag{iCond}==1;
    if(scatterData)
        plot(commonInfo.times{iCond}(inIdx),figureData.data{iCond}(inIdx),'.','color',colors{iCond});
%         if(isfield(figureData,'inOutFlag'))
            %KJ: indicate outliers not used in fit
            %outliers based on value
            outIdx = figureData.inOutFlag{iCond}==0;
            plot(commonInfo.times{iCond}(outIdx),figureData.data{iCond}(outIdx),'o','color',colors{iCond})
            %outliers based on isolation
            outIdx = figureData.inOutFlag{iCond}==-1;
            plot(commonInfo.times{iCond}(outIdx),figureData.data{iCond}(outIdx),'s','color',colors{iCond})
%         end
    end
    
end

lineObj2 = [lineObj{:}];
legend(lineObj2, commonInfo.conditions);
title([figureData.titleBase, ' ', figureData.titleVariable])
xlabel('Time (min)');
ylabel(figureData.yLabel);
ylim([figureData.yMin, figureData.yMax]);

%save and close
if(isempty(outputDirFig2))
    if(~noClose)
        pause;
    end
else
    savefig(figureHandle, [outputDirFig2, filesep, figureData.titleBase, ' ', figureData.titleVariable, '.fig']);
end
if(~noClose)
    close(figureHandle);
end
end

function out = standardErrorFigureError(err,figureData)
if(isempty(err.identifier))
    error('standardErrorFigure:standardErrorFigureError',['Could not plot ' figureData.titleBase ' ' figureData.titleVariable])
end
try
    warning('standardErrorFigure:standardErrorFigureError',['Could not plot ' figureData.titleBase ' ' figureData.titleVariable]);
    disp(err);
    disp(figureData);
    out = [];
    disp('Continuing ...');
catch newError
    disp('standardErrorFigure:ErrorInErrorFunction','An error has occured in timeCourseAnalysis.scatterFigure/scatterFigureError');
    disp(getReport(newError));
end
end

function out = standardErrorFigureErrorParallel(err,d)
% Take advantage of the newer MException class to print more
% informative errors
try
    figureData = d.InputArguments{1};
    warning('standardErrorFigureErrorParallel:standardErrorFigureError',['Could not plot ' figureData.titleBase ' ' figureData.titleVariable]);
    disp(getReport(err));
    disp(d);
    disp(figureData);
    out = [];
    disp('Continuing ...');
catch newError
    warning('standardErrorFigureErrorParallel:ErrorInErrorFunction','An error has occured in timeCourseAnalysis.scatterFigure/scatterFigureErrorParallel');
    try
        disp(getReport(newError));
    catch newNewError
        warning('standardErrorFigureErrorParallel:ErrorInErrorFunction2','A deeper error has occured in timeCourseAnalysis.scatterFigure/scatterFigureErrorParallel');
    end
end
end

