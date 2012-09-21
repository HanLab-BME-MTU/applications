
for i = 1 : length(movieStructAlphaV);
    
    disp(num2str(i))
    
    if movieStructAlphaV(i).activityLevel > 0
        
        tmp = movieStructAlphaV(i).fileName;
        cd([tmp{1} '/analysisAlphaV/furtherAnalysis'])
        
        load tracksDiffusionLength5InMask.mat
        diffModeAnalysisRes = trackDiffModeAnalysis(tracksFinal,diffModeDividerStructIntegrin);
        save('diffusionModeClassification','diffModeAnalysisRes','diffModeDividerStructIntegrin');
        
    end
    
end

