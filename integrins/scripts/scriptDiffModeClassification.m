
for i = 1 : length(movieStructFarnSim20p5);
    
    disp(num2str(i))
    
    if movieStructFarnSim20p5(i).activityLevel > 0
        
        tmp = movieStructFarnSim20p5(i).fileName{1};
        %         cd([tmp '/analysisAlphaV/furtherAnalysis'])
        tmp = regexprep(tmp,'/','\');
        topDir = ['C:\kjData\' tmp(33:end)];
        cd([topDir '\analysisFarnSim20p5\furtherAnalysis'])
        
        load tracksDiffusionLength5InMask.mat
        diffModeAnalysisRes = trackDiffModeAnalysis(tracksFinal,diffModeDividerStructIntegrin);
        save('diffusionModeClassification','diffModeAnalysisRes','diffModeDividerStructIntegrin');
        
    end
    
end

