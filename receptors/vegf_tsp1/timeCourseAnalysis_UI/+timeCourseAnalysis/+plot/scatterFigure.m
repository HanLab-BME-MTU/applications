function [ figureHandles ] = scatterFigure( commonInfo, figureData, outputDirFig, doInParallel )
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

    if(doInParallel)
        figureHandles = pararrayfun_progress( ... 
            @(f) scatterIndividualFigure( ...
                commonInfo, ...
                f, ...
                colors, ...
                outputDirFig, ...
                noClose), ...
            figureData, ...
            'UniformOutput', false, ...
            'ErrorHandler',@scatterFigureError ...
        );
    else
        figureHandles = arrayfun( ... 
            @(f) scatterIndividualFigure( ...
                commonInfo, ...
                f, ...
                colors, ...
                outputDirFig, ...
                noClose), ...
            figureData, ...
            'UniformOutput', false, ...
            'ErrorHandler',@scatterFigureError ...
        );
    end
    
    if(nargout > 0)
        figureHandles = [figureHandles{:}];
    end
        
end
function figureHandle = scatterIndividualFigure(commonInfo, figureData, colors, outputDirFig,noClose)

    plotTitle = [figureData.titleBase ' ' figureData.titleVariable];


    figureHandle = figure('Name', plotTitle);
    hold on;
    %plots all data and stores all line handles
    lineHandle = cellfun(@plot, figureData.fitData, commonInfo.times, figureData.data, 'UniformOutput', false);
    %KJ: indicate outliers not used in fit
    if(~isempty(figureData.inOutFlag))
        outIdx = cellfun(@(inOutFlag) ~inOutFlag,figureData.inOutFlag,'UniformOutput',false);
        cellfun(@(times,data,outIdx) plot(times(outIdx),data(outIdx),'ko'),commonInfo.times,figureData.data,outIdx,'UniformOutput',false);
    end
    %plot vertical lines indicating aligning Times
    nCond = numel(figureData.data);
    for iCond = 1:nCond
        plot([commonInfo.timeShift(iCond), commonInfo.timeShift(iCond)], [figureData.yMax, figureData.yMin], 'Color', colors{iCond});
    end
    %change the color so that color of data and fit match
    cellfun(@(x, y) set(x, 'Color', y), lineHandle, colors);
    %create legends only contain the fit
    fitHandle = [lineHandle{:}];
    legend(fitHandle(2,:), commonInfo.conditions);
    %axis limit
    if ~isempty(figureData.yMax)
        ylim([figureData.yMin, figureData.yMax]);
    end
    %label axis
    xlabel('Time (min)');
    ylabel(figureData.yLabel);
    title(plotTitle);
    %% Saving
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
end
function out = scatterFigureError(err,figureData)
    if(isempty(err.identifier))
        error('scatterFigure:scatterFigureError',['Could not plot ' figureData.titleBase ' ' figureData.titleVariable])
    end
    warning('scatterFigure:scatterFigureError',['Could not plot ' figureData.titleBase ' ' figureData.titleVariable]);
    disp(err);
    disp(figureData);
    out = [];
end