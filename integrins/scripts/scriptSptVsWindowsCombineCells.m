
j = 0;

for i = 1 : length(movieStructAlphaV)
    
    activityLevel = movieStructAlphaV(i).activityLevel;
    
    if activityLevel > 1
        
        j = j + 1;
        
        tmp = movieStructAlphaV(i).fileName{1};
        %         tmp = tmp(11:end);
        cd([tmp '/analysisAlphaV/furtherAnalysis/adaptiveWindows'])
        
        load particleBehaviorAdaptiveWindows121002.mat
        sptPropInWindowAlphaVInd(:,j) = sptPropInWindow;
        windowDistFromEdgeAlphaVInd(:,j) = windowDistFromEdge;
        
    end
    
end

cd /home/kj35/files/LCCB/receptors/Galbraiths/analysis/120907_MovieInfoStructures
