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

%get number of conditions
nCond = length(subData);

%determine maximum y value to determine y axis limit
finiteData = cellfun(@(x) x(isfinite(x(:))),subData,'UniformOutput',false);
finiteData = cellfun(@(x) x(:),finiteData,'UniformOutput',false);
finiteData = vertcat(finiteData{:});
maxValue = max(finiteData);
if ~isempty(maxValue)
    axisTick = 10^(round(log10(maxValue) + 0.2)-1);
else
    axisTick = 10^(round(log10(1) + 0.2)-1);
end
if isYMin0
    yMin = 0;
else
    minValue = min(finiteData);
    if ~isempty(minValue)
        axisTick = max(10^(round(log10(abs(minValue)) + 0.2)-1), axisTick);
        yMin = (ceil(minValue / axisTick) - 1) * axisTick;
    else
        axisTick = max(10^(round(log10(abs(0)) + 0.2)-1), axisTick);
        yMin = (ceil(0 / axisTick) - 1) * axisTick;
    end
end
yMax = (floor(maxValue / axisTick) + 1) * axisTick;

%plot by column
figureData(nColumns) = struct;
for iColumns = 1:nColumns
    
    subSubData = cellfun(@(x) x(:,iColumns), subData, 'UniformOutput', false);
    plotTitle = [title_Base ' ' title_Variable{iColumns}];
    
    %KJ: identify inlier and outlier data points
    inOutFlag = cellfun(@(x) ones(size(x)),subSubData,'UniformOutput',false);
    if(commonInfo.parameters.detectOutliers_k_sigma > 0)
        for iData = 1 : length(subSubData)
            [outIdx,inIdx] = detectOutliers(subSubData{iData},commonInfo.parameters.detectOutliers_k_sigma);
            inOutFlag{iData}(outIdx) = 0; %outliers get flag 0
            inOutFlag{iData}(inIdx)  = 1; %inliers get flag 1
        end
    end
    
    %KJ: average over time intervals
    %designate isolated datapoints (i.e. those in time intervals with very
    %few datapoints) as outliers
    [dataAve, inOutFlag] = calcMultipleDataAve(subSubData, commonInfo.times, inOutFlag, commonInfo.parameters.aveInterval, commonInfo.parameters.shiftTime);
    for iCond = 1 : nCond
        if(all(isnan(dataAve{iCond}.dataAve)))
            warning(['Could not average data for "', plotTitle, '" for CML ', num2str(iCond)]);
        end
    end
    
    %spline fit
    [fitData] = calcMultipleSmoothingSpline(subSubData, commonInfo.times, inOutFlag, commonInfo.parameters.smoothingPara);
    for iCond = 1 : nCond
        if(isempty(fitData{iCond}))
            warning(['Could not calculate smoothing spline for "', plotTitle '" for CML ', num2str(iCond)]);
        end
    end
    
    %Save data in figureData
    figureData(iColumns).titleBase = title_Base;
    figureData(iColumns).titleVariable = title_Variable{iColumns};
    figureData(iColumns).fitData = fitData;
    figureData(iColumns).dataAve = dataAve;
    figureData(iColumns).data = subSubData;
    figureData(iColumns).figureDir = [commonInfo.outputDirFig filesep plotTitle '.fig'];
    figureData(iColumns).yMax = yMax;
    figureData(iColumns).yMin = yMin;
    figureData(iColumns).yLabel = yLabelName;
    figureData(iColumns).inOutFlag = inOutFlag;
    
end

end
