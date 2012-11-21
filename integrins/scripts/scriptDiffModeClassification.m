
for i = 1 : length(movieStructAlphaV724Trunc);
    
    disp(num2str(i))
    
    if movieStructAlphaV724Trunc(i).activityLevel > 0
        
        tmp = movieStructAlphaV724Trunc(i).fileName;
        cd([tmp{1} '/analysisAlphaV724Trunc/furtherAnalysis'])
        
        load tracksDiffusionLength5InMask.mat
        diffModeAnalysisRes = trackDiffModeAnalysis(tracksFinal,diffModeDividerStructIntegrin);
        save('diffusionModeClassification','diffModeAnalysisRes','diffModeDividerStructIntegrin');
        
    end
    
end

