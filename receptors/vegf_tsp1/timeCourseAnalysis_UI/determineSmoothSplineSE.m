function [fitError] = determineSmoothSplineSE(data, time, nBoot, timeResolution, timeLimit, smoothingPara, inOutFlag)
%determines the standard error of a fitted curve using Bootstrapping method
%
%SYNOPSIS [fitError] = determineSE_Bootstrp(data, time, nBoot, timeResolution)
%
%INPUT
%   data            : cellArray of data sets in a column
%   time            : cellArray of time points in a column
%                     The cell element of time corresponds to cell element
%                     of data
%   nBoot           : number of bootstrap data sets to use for bootstrap
%                     analysis.
%   timeResolution  : time resolution of bootstrap analysis
%   timeLimit       : [min, max] limit of analysis
%   inOutFlag       : cellArray of inlier and outlier flags in a column
%
%OUTPUT
%   fitError    : cellArray of standard error

%% Analysis
%Divide data up based on conditions
fitError = cellfun(@doAnalysis, data, time, timeLimit, inOutFlag, 'UniformOutput', false);
%nested function that deals with each data set
    function fitStd = doAnalysis(subData, subTime, subTimeLimit, inOutFlag)
        mask = ~(isnan(subData) | isinf(subData));
        subData = subData(mask);
        subTime = subTime(mask);
        inOutFlag = inOutFlag(mask);
        %KJ: discard outliers not used in original fit
        inIdx = find(inOutFlag == 1);
        subTime = subTime(inIdx);
        subData = subData(inIdx);
        nData = numel(subData);
        %create random numbers (bootstraping)
        newDataIndx = randi(nData, nBoot, nData);
        newDataIndx = sort(newDataIndx, 2);
        mask = true(1,nBoot);
        fitData = cell(1,nBoot);
        parfor iBoot = 1:nBoot
            %Create new data set from the old by picking randomly
            newData = subData(newDataIndx(iBoot, :));
            newTime = subTime(newDataIndx(iBoot, :));
            %smooth out the data
            smoothData = smooth(newData, 5);
            smoothData = smoothData(2:end-2);
            smoothTime = smooth(newTime, 5);
            smoothTime = smoothTime(2:end-2);
            %fit data
            try
                fitData{iBoot} = fit(smoothTime, smoothData, 'smoothingSpline', 'smoothingParam', smoothingPara);
            catch
                mask(iBoot) = false;
            end
        end
        fitData = fitData(mask);
        fitStd = arrayfun(@(x) std(cellfun(@(y) y(x), fitData)), subTimeLimit(1):timeResolution:subTimeLimit(2));
    end
end

