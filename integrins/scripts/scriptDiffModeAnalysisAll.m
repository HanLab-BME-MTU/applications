
tracksAll = [];
for i=1:length(analysisAlphaVTruncInd)
    tmp = analysisAlphaVTruncInd(i).fileName;
    load([tmp{1} '/analysisAlphaVTrunc/furtherAnalysis/tracksDiffusionLength5InMask.mat'])
    tracksAll = [tracksAll; tracksFinal];
end

[modeParamAll,expParamAll] = getDiffModes(tracksAll,5,0.05,1,10,2,'AlphaV All');
