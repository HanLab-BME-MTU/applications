function figureData = calcFigureData(commonInfo, subData, title_Base, title_Variable, yLabelName, isYMin0)
% timeCourseAnalysis.calcFigureData
%
% INPUT
% commonInfo:      structure
%   .parameters:    structure
%       .detectOutliers_k_sigma:  k parameter for detectOutliers
%       .smoothingPara smoothing: parameter for calcMultipleSmoothingSpline
%   .times: times for calcMultipleSmoothingSpline
%   .ouputDirFig: directory to save figure
% subData:         numeric matrix where columns represent data
% title_Base:      common title for all columns
% title_Variable:  column labels for the columns in subData
% yLabelName: label for the yaxis
% isYMin0: set the minimum y-value to zero
%
% OUTPUT
% figureData
%
% Mark Kittisopikul, November 2015
    if(nargin < 5)
        yLabelName = '';
    end
    if(nargin < 6)
        isYMin0 = true;
    end
    %initialization
    nColumns = numel(title_Variable);
    %determine maximum y value to determine y axis limit
    finiteData = cellfun(@(x) x(isfinite(x(:))),subData,'UniformOutput',false);
    finiteData = vertcat(finiteData{:});
    maxValue = max(finiteData);
    axisTick = 10^(round(log10(maxValue) + 0.2)-1);
    if isYMin0
        yMin = 0;
    else
        minValue = min(finiteData);
        axisTick = max(10^(round(log10(abs(minValue)) + 0.2)-1), axisTick);
        yMin = (ceil(minValue / axisTick) - 1) * axisTick;
    end
    yMax = (floor(maxValue / axisTick) + 1) * axisTick;
    %plot by column
    figureData(nColumns) = struct;
    for iColumns = 1:nColumns
        subSubData = cellfun(@(x) x(:,iColumns), subData, 'UniformOutput', false);
        %KJ: identify inlier and outlier data points
        inOutFlag = cellfun(@(x) true(size(x)),subSubData,'UniformOutput',false);
        if(commonInfo.parameters.detectOutliers_k_sigma > 0)
            for iData = 1 : length(subSubData)
                [outIdx,inIdx] = detectOutliers(subSubData{iData},3);
                inOutFlag{iData}(outIdx) = false; %outliers get flag 0
                inOutFlag{iData}(inIdx)  = true; %inliers get flag 1
            end
        end
        plotTitle = [title_Base ' ' title_Variable{iColumns}];
        [fitData] = calcMultipleSmoothingSpline(subSubData, commonInfo.times, inOutFlag, commonInfo.parameters.smoothingPara);
        if(isempty(fitData))
            warning(['Could not calculate smoothing spline for ' plotTitle]);
        end
        %Save data in figureData
        figureData(iColumns).titleBase = title_Base;
        figureData(iColumns).titleVariable = title_Variable{iColumns};
        figureData(iColumns).fitData = fitData;
        figureData(iColumns).data = subSubData;
        figureData(iColumns).figureDir = [commonInfo.outputDirFig filesep plotTitle '.fig'];
        figureData(iColumns).yMax = yMax;
        figureData(iColumns).yMin = yMin;
        figureData(iColumns).yLabel = yLabelName;
        figureData(iColumns).inOutFlag = inOutFlag;
    end
end
