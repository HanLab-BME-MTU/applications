
for i = 1 : length(movieStructAlphaVPax);
    
    disp(num2str(i))
    
    if movieStructAlphaVPax(i).activityLevel > 0
        
        tmp = movieStructAlphaVPax(i).fileName;
        cd([tmp{1} '/analysisAlphaV/furtherAnalysis'])
        
        load tracksDiffusionLength5InMask.mat
        diffModeAnalysisRes = trackDiffModeAnalysis(tracksFinal,diffModeDividerStructIntegrin);
        save('diffusionModeClassification','diffModeAnalysisRes','diffModeDividerStructIntegrin');
        
    end
    
end

