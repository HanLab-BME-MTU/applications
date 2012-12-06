
for i = 1 : length(movieStructAlphaVFixed);
    
    disp(num2str(i))
    
    if movieStructAlphaVFixed(i).activityLevel > 0
        
        tmp = movieStructAlphaVFixed(i).fileName;
        cd([tmp{1} '/analysisAlphaV/furtherAnalysis'])
        
        load tracksDiffusionLength5InMask.mat
        diffModeAnalysisRes = trackDiffModeAnalysis(tracksFinal,diffModeDividerStructIntegrin);
        save('diffusionModeClassification','diffModeAnalysisRes','diffModeDividerStructIntegrin');
        
    end
    
end

