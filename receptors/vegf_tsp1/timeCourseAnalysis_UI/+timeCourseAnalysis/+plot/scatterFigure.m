function [ figureHandles ] = scatterFigure( commonInfo, figureData, outputDirFig )
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
    if(nargout > 0)
        noClose = true;
    else
        noClose = false;
    end
    
    colors = timeCourseAnalysis.plot.getColors(commonInfo.conditions);

    figureHandles = arrayfun( ... 
        @(f) scatterIndividualFigure( ...
            f.fitData, ...
            commonInfo.times, ...
            f.data, ...
            commonInfo.timeShift, ...
            f.yMax, ...
            f.yMin, ...
            colors, ...
            commonInfo.conditions, ...
            f.yLabel, ...
            [f.titleBase ' ' f.titleVariable], ...
            outputDirFig, ...
            noClose), ...
        figureData, ...
        'UniformOutput', false ...
    );
    
    if(nargout > 0)
        figureHandles = [figureHandles{:}];
    end
        
end
function figureHandle = scatterIndividualFigure(dataFit, times, data, alignTimes, yMax, yMin, colors, names, yLabelName, plotTitle, outputDirFig,noClose)
    figureHandle = figure('Name', plotTitle);
    hold on;
    %plots all data and stores all line handles
    lineHandle = cellfun(@plot, dataFit, times, data, 'UniformOutput', false);
    %plot vertical lines indicating aligning Times
    nCond = numel(data);
    for iCond = 1:nCond
        plot([alignTimes(iCond), alignTimes(iCond)], [yMax, yMin], 'Color', colors{iCond});
    end
    %change the color so that color of data and fit match
    cellfun(@(x, y) set(x, 'Color', y), lineHandle, colors);
    %create legends only contain the fit
    fitHandle = [lineHandle{:}];
    legend(fitHandle(2,:), names);
    %axis limit
    if ~isempty(yMax)
        ylim([yMin, yMax]);
    end
    %label axis
    xlabel('Time (min)');
    ylabel(yLabelName);
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
