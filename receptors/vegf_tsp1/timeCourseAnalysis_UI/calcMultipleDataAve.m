function [dataAve, inOutFlag] = calcMultipleDataAve(data, times, inOutFlag, aveInterval, shiftTime)
%Scatter plots given data sets in one plot and averages in aveInterval intervals
%
%SYNOPSIS [dataAve] = calcMultipleDataAve(data, times, inOutFlag, aveInterval)
%
%INPUT
%   data            : cell array of data sets (each set must be a column)
%   times           : cell array of times corresponding to the dataset (each set must be a column)
%   plotTitle       : name of the plot
%   inOutFlag       : cell array of flags indicating inlier data points
%                     (value 1) and outlier data points (value 0)
%   aveInterval     : Interval for averaging. Default: 3 min
%   shiftTime       : Value by which time is shifted. Needed to shift time
%                     back temporarily for this analysis.
%
%OUTPUT
%   aveData         :
%
% Khuloud Jaqaman, March 2017

%% Initialization
%assign default value
if isempty(aveInterval)
    aveInterval = 3;
end
if isempty(shiftTime)
    shiftTime = zeros(size(data));
end

shiftTimeCell = cell(size(data));
for iCell = 1 : length(shiftTimeCell)
    shiftTimeCell{iCell} = shiftTime(iCell);
end

%creates figure and stores the figure handle
fxn = @(varargin) calcMultipleDataAvePerCondition(varargin{:}, aveInterval);
[dataAve, inOutFlag] = cellfun(fxn, data, times, inOutFlag, shiftTimeCell, 'UniformOutput',false);

end


%% sub-function

function [dataAve, inOutFlag] = calcMultipleDataAvePerCondition(data, times, inOutFlag, shiftTime, aveInterval)

try
    %     %remove nan and inf
    %     mask = isfinite(data);
    %     data = data(mask);
    %     times = times(mask);
    %     inOutFlag = inOutFlag(mask);
    %
    %     %KJ: discard outliers
    %     timesIn = times(inOutFlag);
    %     dataIn = data(inOutFlag);
    
    %divide time points into bins and take time and data average and
    %std
    binID = round((times-shiftTime)/aveInterval);
    binUnique = unique(binID);
    nBin = length(binUnique);
    [nDataPoints,dataInAve,dataInStd,timeInAve,timeInStd] = deal(NaN(nBin,1));
    for iBin = 1 : nBin
        indx = find(binID == binUnique(iBin));
        finiteFlagBin = isfinite(data(indx));
        inOutFlagBin = inOutFlag(indx);
        indxGood = indx(finiteFlagBin&inOutFlagBin==1);
        nDP = length(indxGood);
        nDataPoints(iBin) = nDP;
        if nDP >= 5
            dataInAve(iBin) = mean(data(indxGood));
            timeInAve(iBin) = mean(times(indxGood));
            dataInStd(iBin) = std(data(indxGood));
            timeInStd(iBin) = std(times(indxGood));
        else
            inOutFlag(indx) = -1;
        end
    end
    
    %output
    dataAve = struct('timeAve',timeInAve,'dataAve',dataInAve,...
        'timeStd',timeInStd,'dataStd',dataInStd,'nDataPts',nDataPoints);
    
catch
    dataAve = [];
end

end