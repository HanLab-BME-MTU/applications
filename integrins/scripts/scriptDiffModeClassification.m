
for i = 1 : length(movieStructAlphaVY773A);
    
    disp(num2str(i))
    
    if movieStructAlphaVY773A(i).activityLevel > 0
        
        tmp = movieStructAlphaVY773A(i).fileName{1};
        %         cd([tmp '/analysisAlphaVY773A/furtherAnalysis'])
        tmp = regexprep(tmp,'/','\');
        topDir = ['C:\kjData\' tmp(33:end)];
        cd([topDir '\analysisAlphaVY773A\furtherAnalysis'])
        
        load tracksDiffusionLength5InMask.mat
        diffModeAnalysisRes = trackDiffModeAnalysis(tracksFinal,diffModeDividerStructIntegrin);
        save('diffusionModeClassification','diffModeAnalysisRes','diffModeDividerStructIntegrin');
        
    end
    
end
