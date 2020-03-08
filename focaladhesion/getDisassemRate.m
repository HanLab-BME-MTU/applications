function [disassemRate,bestModel,bestSummary] = getDisassemRate(tRangeMin,TS,minLength)
%function assemRate = getAssemRate(tRangeMin,curAmpTotal) calculates
%assembly rate using Webb 2004 NCB method.
% input:
%       tRangeMin           time series of time in minute
%       TS                  time series of signal
% output:
%       assemRate           assembly rate
%       bestModel           The best linear model out of fitting
% Sangyon Han March 4, 2020

%% Get the maximum amp, time range
if nargin<3
    minLength=9;
end

splineParam=0.01; 
sd_spline= csaps(tRangeMin,TS,splineParam);
sd=ppval(sd_spline,tRangeMin);
% Find the maximum
[~,maxSdInd] = max(sd);
% maxAmpFrame = tRangeMin(maxSdInd);
% minLength = 9;
maxSdInd = min(maxSdInd,length(TS)-minLength+1);

% nSampleStart=min(minLength,floor((maxSdInd)/2));

% minAmp = min(curAmpTotal);
% maxAmp = max(curAmpTotal);

pp=0;
fitSummary(length(TS)-minLength-maxSdInd+2)=struct('adjRsquared',NaN,'rSquared',NaN,'pValue',NaN,'slope',NaN,'time_length',NaN,'image_count',NaN,'startingTS',NaN); 
statModel = cell(1,length(TS)-minLength-maxSdInd+2);
for ii= minLength:length(TS)-maxSdInd
    pp=pp+1;
    curTS = TS(end-ii:end);
    curTSnorm = curTS(1)./curTS;
    curTSnorm(curTSnorm < 0) = NaN;
    curTRange = tRangeMin(end-ii:end);

    statModel{pp} = fitlm(curTRange, log(curTSnorm));
    % thisModel2 = fitlm(tRangeMin(1:ii), TS(1:ii));

    fitSummary(pp).adjRsquared = statModel{pp}.Rsquared.Adjusted;
    fitSummary(pp).rSquared = statModel{pp}.Rsquared.Ordinary;
    fitSummary(pp).pValue = statModel{pp}.Coefficients.pValue;
    fitSummary(pp).slope = statModel{pp}.Coefficients.Estimate(2);
    fitSummary(pp).time_length = tRangeMin(ii);
    fitSummary(pp).image_count = ii;
    fitSummary(pp).startingTS = curTS(1);
end

%% Find the best models
adjRS_all = arrayfun(@(x) x.adjRsquared,fitSummary);
[~,maxRframe] = max(adjRS_all);
bestSummary = fitSummary(maxRframe);
bestModel = statModel{maxRframe};
disassemRate = bestSummary.slope;



% nSampleEnd=min(9,floor((length(tRange)-maxSdInd)*2/3));
% sigTtest=ttest2(curAmpTotal(end-nSampleEnd:end),curAmpTotal(maxSdInd:maxSdInd+nSampleEnd));
% if ~isnan(sigTtest)
%     if nSampleEnd>4 && sigTtest && ...
%             mean(curAmpTotal(end-nSampleEnd:end))<mean(curAmpTotal(maxSdInd:maxSdInd+nSampleEnd))
%         [~,disassemRate] = regression(tIntervalMin*(tRange(maxSdInd:end)),...
%             log(curAmpTotal(maxSdInd) ./curAmpTotal(maxSdInd:end)));
%     else
%         disassemRate = NaN;
%     end
% else
%     disassemRate = NaN;
% end



end

