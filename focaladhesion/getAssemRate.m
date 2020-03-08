function [assemRate,bestModel,bestSummary] = getAssemRate(tRangeMin,TS,minLength)
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
splineParam=0.01; 
sd_spline= csaps(tRangeMin,TS,splineParam);
sd=ppval(sd_spline,tRangeMin);
% Find the maximum
[~,maxSdInd] = max(sd);
% maxAmpFrame = tRangeMin(maxSdInd);
% minLength = 9;
maxSdInd=max(maxSdInd,minLength+1);

% nSampleStart=min(minLength,floor((maxSdInd)/2));

% minAmp = min(curAmpTotal);
% maxAmp = max(curAmpTotal);

TSnorm = TS/TS(1);
TSnorm(TSnorm < 0) = NaN;

pp=0;
fitSummary(maxSdInd-minLength+1)=struct('adjRsquared',NaN,'rSquared',NaN,'pValue',NaN,'slope',NaN,'time_length',NaN,'image_count',NaN); 
statModel = cell(1,maxSdInd-minLength+1);
for ii= minLength:maxSdInd
    pp=pp+1;
    statModel{pp} = fitlm(tRangeMin(1:ii), log(TSnorm(1:ii)));
    % thisModel2 = fitlm(tRangeMin(1:ii), TS(1:ii));

    fitSummary(pp).adjRsquared = statModel{pp}.Rsquared.Adjusted;
    fitSummary(pp).rSquared = statModel{pp}.Rsquared.Ordinary;
    fitSummary(pp).pValue = statModel{pp}.Coefficients.pValue;
    fitSummary(pp).slope = statModel{pp}.Coefficients.Estimate(2);
    fitSummary(pp).time_length = tRangeMin(ii);
    fitSummary(pp).image_count = ii;
end

%% Find the best models
adjRS_all = arrayfun(@(x) x.adjRsquared,fitSummary);
[~,maxRframe] = max(adjRS_all);
bestSummary = fitSummary(maxRframe);
bestModel = statModel{maxRframe};
assemRate = bestSummary.slope;


    
    
% if ~all(isnan(TS(maxSdInd-nSampleStart+1:maxSdInd))) && ~all(isnan(TS(1:nSampleStart)))
%     sigTtest = ttest2(TS(1:nSampleStart),TS(maxSdInd-nSampleStart+1:maxSdInd));
%     if ~isnan(sigTtest)
%         if nSampleStart>4 && sigTtest && ...
%                 nanmean(TS(1:nSampleStart))<nanmean(TS(maxSdInd-nSampleStart+1:maxSdInd))
%             [~,assemRate] = regression(tIntervalMin*(tRange(1:maxSdInd)),...
%                 log(curTrack.ampTotal(curTrack.startingFrameExtra:maxAmpFrame)/...
%                 curTrack.ampTotal(curTrack.startingFrameExtra)));
%         else
%             assemRate = NaN;
%         end
%     else
%         assemRate = NaN;
%     end
% else
%     assemRate = NaN;
% end



end

