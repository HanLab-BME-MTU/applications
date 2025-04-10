function [assemRate,bestModel,bestSummary,tRangeSelected] = getAssemRate(tRangeMin,TS,minLength)
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
N = maxSdInd-minLength+1;
Ncumul = N*(N+1)/2;
fitSummary(Ncumul)=struct('adjRsquared',NaN,'rSquared',NaN,'pValue',NaN,'slope',NaN,'time_length',NaN,'image_count',NaN); 
statModel = cell(1,Ncumul);
for ii= minLength:maxSdInd %This was changing only the last point, not the start point. 
    for jj=1:(ii-minLength+1)
        pp=pp+1;
        TSnorm = TS/TS(jj);
        statModel{pp} = fitlm(tRangeMin(jj:ii), log(TSnorm(jj:ii)));
        % thisModel2 = fitlm(tRangeMin(1:ii), TS(1:ii));

        fitSummary(pp).adjRsquared = statModel{pp}.Rsquared.Adjusted;
        fitSummary(pp).rSquared = statModel{pp}.Rsquared.Ordinary;
        fitSummary(pp).pValue = statModel{pp}.Coefficients.pValue(2);
        fitSummary(pp).slope = statModel{pp}.Coefficients.Estimate(2);
        fitSummary(pp).time_length = tRangeMin(ii);
        fitSummary(pp).image_count = ii;
    end
end

%% Find the best models
adjRS_all = arrayfun(@(x) x.adjRsquared,fitSummary);
[~,maxRframe] = max(adjRS_all);
p_all = arrayfun(@(x) x.pValue,fitSummary);
[~,minPframe] = min(p_all);
if minPframe > maxRframe % We choose frame that is longer even if slope is lower. 
    % This way we prevent to capture insanely high rate
    chosenFrame = minPframe;
else
    chosenFrame = maxRframe;
end
bestSummary = fitSummary(chosenFrame);
bestModel = statModel{chosenFrame};
assemRate = bestSummary.slope;

if assemRate<0
    bestSummary = [];
    bestModel = [];
    assemRate = NaN;
    tRangeSelected=[];
    return
end

pp=0;
for ii= minLength:maxSdInd %This was changing only the last point, not the start point. 
    for jj=1:(ii-minLength+1)
        pp=pp+1;
        if pp==chosenFrame
            startFrame = jj;
            endFrame = ii;
            break
        end
    end
end
% iiRange = minLength:maxSdInd;
tRangeSelected=startFrame:endFrame;

    
    
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

