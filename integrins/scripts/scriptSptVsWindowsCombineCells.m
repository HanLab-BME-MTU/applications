
j = 0;

for i = 1 : length(movieStructFarn)
    
    activityLevel = movieStructFarn(i).activityLevel;
    
    if activityLevel > 1
        
        j = j + 1;
        
        tmp = movieStructFarn(i).fileName{1};
        %         tmp = tmp(11:end);
        cd([tmp '/analysisFarn/furtherAnalysis/adaptiveWindows'])
        
        load particleBehaviorAdaptiveWindows121112.mat
        sptPropInWindowFarnInd(:,j) = sptPropInWindow;
        windowDistFromEdgeFarnInd(:,j) = windowDistFromEdge;
        
    end
    
end

cd /home/kj35/files/LCCB/receptors/Galbraiths/analysis/121203_sptVsWindowsCombCellsNoJumps