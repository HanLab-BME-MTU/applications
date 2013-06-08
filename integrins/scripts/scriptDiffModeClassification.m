
for i = 1 : length(movieStructBeta3);
    
    disp(num2str(i))
    
    if movieStructBeta3(i).activityLevel > 0
        
        tmp = movieStructBeta3(i).fileName{1};
        %         cd([tmp '/analysisBeta3/furtherAnalysis'])
        tmp = regexprep(tmp,'/','\');
        topDir = ['C:\kjData\' tmp(33:end)];
        cd([topDir '\analysisBeta3\furtherAnalysis'])
        
        load tracksDiffusionLength5InMask.mat
        diffModeAnalysisRes = trackDiffModeAnalysis(tracksFinal,diffModeDividerStructIntegrin);
        save('diffusionModeClassification','diffModeAnalysisRes','diffModeDividerStructIntegrin');
        
    end
    
end
