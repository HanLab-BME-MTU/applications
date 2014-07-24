
for i = 1 : length(movieStructLifeact);
    
    disp(num2str(i))
    
    if movieStructLifeact(i).activityLevel > 0
        
        tmp = movieStructLifeact(i).fileName{1};
        tmp = regexprep(tmp,'/','\');
        topDir = ['C:\kjData\' tmp(33:end)];
        cd([topDir '\analysisLifeact\furtherAnalysis'])
        
        load tracksDiffusionLength5InMask.mat
        load diffusionModeClassification.mat
        
        [compMode2MSS,compModeWithMSS2ModeNoMSS] = compDiffAnalysisModeMSS(...
            tracksFinal,diffAnalysisRes,diffModeAnalysisRes,[5 99]);

        save('compDiffAnalysis','compMode2MSS','compModeWithMSS2ModeNoMSS');
        compModeMSSAll(i) = compMode2MSS;
        compModeWithMSS2ModeNoMSSAll(i) = compModeWithMSS2ModeNoMSS;

    end
    
end
