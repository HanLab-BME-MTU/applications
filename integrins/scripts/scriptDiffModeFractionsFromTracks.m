
numMode = 2;
fracModeLifeact = NaN(numMode,length(analysisLifeactInd));

for i = 1 : length(analysisLifeactInd)
    
    tmp = analysisLifeactInd(i).fileName;
    cd([tmp{1} '/analysisLifeact/furtherAnalysis'])
    
    load diffusionModeClassification
    
    diffMode = vertcat(diffModeAnalysisRes.diffMode);
    
    n = hist(diffMode,1:numMode);
    
    fracModeLifeact(:,i) = n / sum(n);
    
    
end

