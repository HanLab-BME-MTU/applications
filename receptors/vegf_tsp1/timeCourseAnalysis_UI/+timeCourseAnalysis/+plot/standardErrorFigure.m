function [ figObj ] = standardErrorFigure( commonInfo, figureData, scatterData, outputDirFig2 )
%standardErrorFigure Plots time course analysis figure with standard error
%
% INPUT
% commonInfo
% figureData
% scatterData - logical whether to scatter data (default: true)
% outputDirFig2 - directory to save the figure
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

    if(nargout > 0)
        noClose = true;
    else
        noClose = false;
    end

nConditions = numel(commonInfo.conditions);
nFig = numel(figureData);
lineObj = cell(1,nConditions);
colors = timeCourseAnalysis.plot.getColors(commonInfo.conditions);

figObj(nFig) = matlab.graphics.GraphicsPlaceholder;

for iFig = 1:nFig
    %open up the figure and make it current figure
    %figObj = openfig(figureData(iFig).figureDir);
    figObj(iFig) = figure();
    hold on;
    %plot standard errors
    for iCond = 1:nConditions
        fitValues = figureData(iFig).fitData{iCond}(commonInfo.analysisTimes{iCond});
        fitValues = fitValues';
        plot(commonInfo.analysisTimes{iCond}, fitValues + figureData(iFig).fitError{iCond}, ':', 'color', colors{iCond});
        plot(commonInfo.analysisTimes{iCond}, fitValues - figureData(iFig).fitError{iCond}, ':', 'color', colors{iCond});
        lineObj{iCond} = plot(commonInfo.analysisTimes{iCond}, fitValues, 'color', colors{iCond});
        plot([commonInfo.timeShift(iCond), commonInfo.timeShift(iCond)], [figureData(iFig).yMax, figureData(iFig).yMin], 'Color', colors{iCond});
        if(scatterData)
            plot(commonInfo.times{iCond},figureData(iFig).data{iCond},'.','color',colors{iCond});
        end
    end
    lineObj2 = [lineObj{:}];
    legend(lineObj2, commonInfo.conditions);
    title([figureData(iFig).titleBase, ' ', figureData(iFig).titleVariable])
    xlabel('Time (min)');
    ylabel(figureData(iFig).yLabel);
    ylim([figureData(iFig).yMin, figureData(iFig).yMax]);
    %save and close
    if(nargin < 4 || isempty(outputDirFig2))
        if(~noClose)
            pause;
        end
    else
        savefig(figObj(iFig), [outputDirFig2, filesep, figureData(iFig).titleBase, ' ', figureData(iFig).titleVariable, '.fig']);
    end
    if(nargout == 0)
        if(~noClose)
            close(figObj(iFig));
        end
    end
end

end

