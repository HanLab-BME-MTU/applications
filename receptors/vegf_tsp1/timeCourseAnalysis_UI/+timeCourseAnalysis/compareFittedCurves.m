function [ fitCompare, compareTimes ] = compareFittedCurves( commonInfo, figureData )
%compareFittedCurves Summary of this function goes here
%   Detailed explanation goes here

%% Compare Fitted Curves
%progress display
progressTextMultiple('Comparing fittted curves', numel(figureData));
%pairs up curve indices for comparison----------------
%And
%determines the commonInfo.compareTimes
%which are simply intersection of two commonInfo.analysisTimes sets
%And
%determines the indx of these overlap------------
nConditions = numel(commonInfo.conditions);
nPair = nConditions * (nConditions - 1) / 2;
timeLimitIndx = commonInfo.timeLimitIndx;
curveIndxPair = cell(1, nPair);
commonIndx1 = cell(1, nPair);
commonIndx2 = cell(1, nPair);
compareTimes = cell(1, nPair);
iPair_ = 0;
for indx1 = 1:nConditions-1
    for indx2 = indx1+1:nConditions
        iPair_ = iPair_ + 1;
        curveIndxPair{iPair_} = [indx1, indx2];
        commonMinIndx = max(timeLimitIndx{indx1}(1), timeLimitIndx{indx2}(1));
        commonMaxIndx = min(timeLimitIndx{indx1}(2), timeLimitIndx{indx2}(2));
        %for indx1
        commonIndx1{iPair_} = [commonMinIndx, commonMaxIndx] - timeLimitIndx{indx1}(1) + 1;
        %for indx2
        commonIndx2{iPair_} = [commonMinIndx, commonMaxIndx] - timeLimitIndx{indx2}(1) + 1;
        %indx to value conversion
        compareTimes{iPair_} = commonInfo.analysisTimes{indx1}(commonIndx1{iPair_}(1):commonIndx1{iPair_}(2));
    end
end
%determine p-values at each timepoints = commonInfo.compareTimes
%deals out figureData to callGetP
fitCompare = arrayfun(@(fd) callGetP(fd,curveIndxPair,nConditions,nPair,compareTimes,commonIndx1, commonIndx2), ...
    figureData, 'UniformOutput', false);
% [figureData.fitCompare] = fitCompare{:};


end

%% Accessory Function: deals out pair of curves to getP
function [fitCompare]  = callGetP(FD,curveIndxPair,nConditions,nPair,compareTimes,commonIndx1, commonIndx2)
    if isempty(curveIndxPair)
        pValue = [];
    else
        pValue = cellfun(@(x, y, index1, index2) getPValues(FD.fitData{x(1)}, FD.fitData{x(2)}, ...
                FD.fitError{x(1)}(index1(1):index1(2)), ...
                FD.fitError{x(2)}(index2(1):index2(2)), ...
                y), ...
            curveIndxPair, ...
            compareTimes, ...
            commonIndx1, ...
            commonIndx2, ...
            'UniformOutput', false, ...
            'ErrorHandler', @getPEH);
    end
    fitCompare(nConditions, nConditions) = struct('geoMeanP', [], 'p', [], 'timeIndx', []);
    for iPair = 1:nPair
        fitCompare(curveIndxPair{iPair}(1), curveIndxPair{iPair}(2)).p = pValue{iPair};
        fitCompare(curveIndxPair{iPair}(1), curveIndxPair{iPair}(2)).geoMeanP = geomean(pValue{iPair});
        fitCompare(curveIndxPair{iPair}(1), curveIndxPair{iPair}(2)).timeIndx = iPair;
    end
    progressTextMultiple();
end
%determines p-values
function [pValue] = getPValues(fitData1, fitData2, fitError1, fitError2, times)
    pValue = arrayfun(@(x, y, t) getP(fitData1(t), fitData2(t), x, y), fitError1, fitError2, times);
end
%determines p-value using z-test
function [pValue] = getP(mean1, mean2, SE1, SE2)
    z = abs(mean1 - mean2) ./ sqrt(SE1.^2 + SE2.^2);
    pValue = 1 - diff(normcdf([-z,z]));
end
%error handle for p value determination
function [] = getPEH(varargin)
    warning(['Curve Comparison for figure ' num2str(varargin{1}.index) ' has failed']);
    error(varargin{1});
end

