
for i = 1 : length(movieStructAlphaVPax)
    
    activityLevel = movieStructAlphaVPax(i).activityLevel;
    
    if activityLevel > 0
        
        tmp = movieStructAlphaVPax(i).fileName{1};
        %         tmp = tmp(11:end);
        
        %         cd([tmp '/analysisAlphaVPax/furtherAnalysis/adaptiveWindows/figuresSptVsWindows121112'])
        %
        %         copyfile('*.fig','/home/kj35/files/LCCB/receptors/Galbraiths/analysis/121113_sptVsWindowsIndCells/');
        
        %         cd([tmp '/analysisCellEdgeSmall'])
        %
        %         copyfile('AlphaVPax*.mov','/home/kj35/files/LCCB/receptors/Galbraiths/analysis/121113_masksIndCells/');
        
        cd([tmp '/analysisAlphaV/furtherAnalysis/spatialMap'])
        
        copyfile('*.fig','/home/kj35/files/LCCB/receptors/Galbraiths/analysis/121128_spatialMapsFinalMasks/');

    end
    
end

cd /home/kj35/files/LCCB/receptors/Galbraiths/analysis/121128_spatialMapsFinalMasks
