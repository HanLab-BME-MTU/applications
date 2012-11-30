
j = 0;

for i = 1 : length(movieStructAlphaV724Trunc)
    
    activityLevel = movieStructAlphaV724Trunc(i).activityLevel;
    
    if activityLevel > 1
        
        j = j + 1;
        
        tmp = movieStructAlphaV724Trunc(i).fileName{1};
        %         tmp = tmp(11:end);
        cd([tmp '/analysisAlphaV724Trunc/furtherAnalysis/adaptiveWindows'])
        
        load particleBehaviorAdaptiveWindows121112.mat
        sptPropInWindowAlphaV724TruncInd(:,j) = sptPropInWindow;
        windowDistFromEdgeAlphaV724TruncInd(:,j) = windowDistFromEdge;
        
    end
    
end

cd /home/kj35/files/LCCB/receptors/Galbraiths/analysis/120907_MovieInfoStructures