
numMode = 2;
fracModeLifeactFixed = NaN(length(movieStructLifeactFixed),numMode);

for i = 1 : length(movieStructLifeactFixed)
    
    activityLevel = movieStructLifeactFixed(i).activityLevel;
    
    if activityLevel > 0
        
        tmp = movieStructLifeactFixed(i).fileName;
        cd([tmp{1} '/analysisLifeact/furtherAnalysis'])
        
        load diffusionModeClassification
        
        diffMode = vertcat(diffModeAnalysisRes.diffMode);
        
        n = hist(diffMode,1:numMode);
        
        fracModeLifeactFixed(i,:) = n / sum(n);
        
    end
    
    
end

