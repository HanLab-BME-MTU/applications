
j = 0;

for i = 1 : length(movieStructTalin)
    
    activityLevel = movieStructTalin(i).activityLevel;
    
    if activityLevel > 1
        
        j = j + 1;
        
        tmp = movieStructTalin(i).fileName{1};
        %         tmp = tmp(11:end);
        cd([tmp '/analysisTalin/furtherAnalysis/adaptiveWindows'])
        
        load particleBehaviorAdaptiveWindows120907.mat
        sptPropInWindowIndTalin(:,j) = sptPropInWindow;
        windowDistFromEdgeIndTalin(:,j) = windowDistFromEdge;
        
    end
    
end

cd /home/kj35/files/LCCB/receptors/Galbraiths/analysis/120907_MovieInfoStructures
