
for i = 1 : length(analysisTalinInd);
    
    disp(num2str(i))
    
    tmp = analysisTalinInd(i).fileName;
    cd([tmp{1} '/analysisTalin/furtherAnalysis'])
    
    load tracksDiffusionLength5InMask.mat
    diffModeAnalysisRes = trackDiffModeAnalysis(tracksFinal,diffModeDividerStructIntegrin);
    save('diffusionModeClassification','diffModeAnalysisRes','diffModeDividerStructIntegrin');
    
end

