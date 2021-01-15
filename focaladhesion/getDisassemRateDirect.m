function [disassemRate,bestModel,bestSummary] = getDisassemRateDirect(tRangeMin,TS,tRangeSelected)
%function assemRate = getAssemRate(tRangeMin,curAmpTotal) calculates
%assembly rate using Webb 2004 NCB method but without iteration.
% input:
%       tRangeMin           time series of time in minute
%       TS                  time series of signal
% output:
%       assemRate           assembly rate
%       bestModel           The best linear model out of fitting
%       tRangeSelected      frame range selected
% Sangyon Han March 4, 2020

%% Get the maximum amp, time range
if isempty(tRangeSelected)
    disassemRate=NaN;
    bestModel = [];
    bestSummary.adjRsquared = NaN;
    bestSummary.rSquared = NaN;
    bestSummary.pValue = NaN;
    bestSummary.slope = NaN;
    bestSummary.time_length = NaN;
    bestSummary.image_count = NaN;
    bestSummary.startingTS = NaN;
    return
end
    
curTS = TS(tRangeSelected);
curTSnorm = curTS(1)./curTS;
curTSnorm(curTSnorm < 0) = NaN;
curTRange = tRangeMin(tRangeSelected);

statModel = fitlm(curTRange, log(curTSnorm));
% thisModel2 = fitlm(tRangeMin(1:ii), TS(1:ii));

fitSummary.adjRsquared = statModel.Rsquared.Adjusted;
fitSummary.rSquared = statModel.Rsquared.Ordinary;
fitSummary.pValue = statModel.Coefficients.pValue;
fitSummary.slope = statModel.Coefficients.Estimate(2);
fitSummary.time_length = tRangeMin(tRangeSelected(1));
fitSummary.image_count = 1;
fitSummary.startingTS = curTS(1);

%% Find the best models
bestSummary = fitSummary;
bestModel = statModel;
disassemRate = bestSummary.slope;

end

