
j = 0;

for i = 1 : length(movieStructAlphaV717Trunc)
    
    activityLevel = movieStructAlphaV717Trunc(i).activityLevel;
    
    if activityLevel > 1
        
        j = j + 1;
        
        tmp = movieStructAlphaV717Trunc(i).fileName{1};
        %         tmp = tmp(11:end);
        cd([tmp '/analysisAlphaV717Trunc/furtherAnalysis/adaptiveWindows'])
        
        load particleBehaviorAdaptiveWindows121112.mat
        sptPropInWindowAlphaV717TruncInd(:,j) = sptPropInWindow;
        windowDistFromEdgeAlphaV717TruncInd(:,j) = windowDistFromEdge;
        
    end
    
end

cd /home/kj35/files/LCCB/receptors/Galbraiths/analysis/121113_sptVsWindowsCombCells
