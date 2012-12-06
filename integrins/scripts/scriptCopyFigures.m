
for i = 1 : length(movieStructAlphaVFixed)
    
    activityLevel = movieStructAlphaVFixed(i).activityLevel;
    
    if activityLevel > 0
        
        tmp = movieStructAlphaVFixed(i).fileName{1};
        %         tmp = tmp(11:end);
        
        %         cd([tmp '/analysisAlphaVFixed/furtherAnalysis/adaptiveWindows/figuresSptVsWindows121112'])
        %
        %         copyfile('*.fig','/home/kj35/files/LCCB/receptors/Galbraiths/analysis/121113_sptVsWindowsIndCells/');
        
        %         cd([tmp '/analysisCellEdgeSmall'])
        %
        %         copyfile('AlphaVFixed*.mov','/home/kj35/files/LCCB/receptors/Galbraiths/analysis/121113_masksIndCells/');
        
        cd([tmp '/analysisAlphaV/furtherAnalysis/spatialMap'])
        
        copyfile('*.fig','/home/kj35/files/LCCB/receptors/Galbraiths/analysis/121128_spatialMapsFinalMasks/');

    end
    
end

cd /home/kj35/files/LCCB/receptors/Galbraiths/analysis/121128_spatialMapsFinalMasks
