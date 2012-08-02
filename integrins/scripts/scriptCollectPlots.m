% clear all
% close all
% cd /home/kj35/files/LCCB/receptors/Galbraiths/data/farnesylAndCellEdge/110829_Cs2C1_CHO_Farn/analysisFarn/furtherAnalysis/adaptiveWindows/figuresSptVsWindows
% h(1) = open('fractionMode1_protType1.fig');
% cd /home/kj35/files/LCCB/receptors/Galbraiths/data/farnesylAndCellEdge/110829_Cs1C2_CHO_Farn/analysisFarn/furtherAnalysis/adaptiveWindows/figuresSptVsWindows
% h(2) = open('fractionMode1_protType1.fig');
% cd /home/kj35/files/LCCB/receptors/Galbraiths/data/farnesylAndCellEdge/110829_Cs1C1_CHO_Farn/analysisFarn/furtherAnalysis/adaptiveWindows/figuresSptVsWindows
% h(3) = open('fractionMode1_protType1.fig');
% cd /home/kj35/files/LCCB/receptors/Galbraiths/data/alphaVandCellEdge/110909_Cs2C1_CHO_AV/analysisAlphaV/furtherAnalysis/adaptiveWindows/figuresSptVsWindows/
% h(4) = open('fractionMode1_protType1.fig');
% cd /home/kj35/files/LCCB/receptors/Galbraiths/data/alphaVandCellEdge/110120/Cs1_CHO03/Cs1_CHO03A/analysisAlphaV/furtherAnalysis/adaptiveWindows/figuresSptVsWindows/
% h(5) = open('fractionMode1_protType1.fig');
% cd /home/kj35/files/LCCB/receptors/Galbraiths/data/alphaVandCellEdge/110120/Cs1_CHO02/Cs1_CHO02B/analysisAlphaV/furtherAnalysis/adaptiveWindows/figuresSptVsWindows/
% h(6) = open('fractionMode1_protType1.fig');
% cd /home/kj35/files/LCCB/receptors/Galbraiths/data/alphaVandCellEdge/110120/Cs1_CHO02/Cs1_CHO02A/analysisAlphaV/furtherAnalysis/adaptiveWindows/figuresSptVsWindows/
% h(7) = open('fractionMode1_protType1.fig');
% cd /home/kj35/files/LCCB/receptors/Galbraiths/data/alphaVandCellEdge/110120/Cs1_CHO01/Cs1_CHO01A/analysisAlphaV/furtherAnalysis/adaptiveWindows/figuresSptVsWindows/
% h(8) = open('fractionMode1_protType1.fig');
% cd /home/kj35/files/LCCB/receptors/Galbraiths/data/alphaVandCellEdge/110114/Cs3_CHO03/Cs3_CHO03B/analysisAlphaV/furtherAnalysis/adaptiveWindows/figuresSptVsWindows/
% h(9) = open('fractionMode1_protType1.fig');
% cd /home/kj35/files/LCCB/receptors/Galbraiths/data/alphaVandCellEdge/110114/Cs3_CHO03/Cs3_CHO03A/analysisAlphaV/furtherAnalysis/adaptiveWindows/figuresSptVsWindows/
% h(10) = open('fractionMode1_protType1.fig');
% cd /home/kj35/files/LCCB/receptors/Galbraiths/data/alphaVandCellEdge/110114/Cs3_CHO01/Cs3_CHO01B/analysisAlphaV/furtherAnalysis/adaptiveWindows/figuresSptVsWindows/
% h(11) = open('fractionMode1_protType1.fig');
% cd /home/kj35/files/LCCB/receptors/Galbraiths/data/alphaVandCellEdge/110114/Cs3_CHO01/Cs3_CHO01A/analysisAlphaV/furtherAnalysis/adaptiveWindows/figuresSptVsWindows/
% h(12) = open('fractionMode1_protType1.fig');
% cd /home/kj35/files/LCCB/receptors/Galbraiths/data/alphaVandCellEdge/110114/Cs1_CHO03/Cs1_CHO03B/analysisAlphaV/furtherAnalysis/adaptiveWindows/figuresSptVsWindows/
% h(13) = open('fractionMode1_protType1.fig');
% cd /home/kj35/files/LCCB/receptors/Galbraiths/data/alphaVandCellEdge/110114/Cs1_CHO03/Cs1_CHO03A/analysisAlphaV/furtherAnalysis/adaptiveWindows/figuresSptVsWindows/
% h(14) = open('fractionMode1_protType1.fig');
% cd /home/kj35/files/LCCB/receptors/Galbraiths/data/alphaVandCellEdge/110114/Cs1_CHO02/Cs1_CHO02A/analysisAlphaV/furtherAnalysis/adaptiveWindows/figuresSptVsWindows/
% h(15) = open('fractionMode1_protType1.fig');
% 
% for i=1:15
%     ha(i) = gca(h(i));
% end
% 
% for i=1:15
%     axis(ha(i),[-3.2 20 0 30]);
% end

clear all
close all
cd /home/kj35/files/LCCB/receptors/Galbraiths/data/alphaV717TruncAndCellEdge/120608_Cs3C2_Av717Trunc/analysisAlphaV717/furtherAnalysis/adaptiveWindows/figuresSptVsWindows
h(1) = open('fractionMode1_protType1.fig');
cd /home/kj35/files/LCCB/receptors/Galbraiths/data/alphaV717TruncAndCellEdge/120608_Cs3C1_Av717Trunc/analysisAlphaV717/furtherAnalysis/adaptiveWindows/figuresSptVsWindows
h(2) = open('fractionMode1_protType1.fig');
cd /home/kj35/files/LCCB/receptors/Galbraiths/data/alphaV717TruncAndCellEdge/120608_Cs2C2_Av717Trunc/analysisAlphaV717/furtherAnalysis/adaptiveWindows/figuresSptVsWindows
h(3) = open('fractionMode1_protType1.fig');
cd /home/kj35/files/LCCB/receptors/Galbraiths/data/alphaV717TruncAndCellEdge/120608_Cs1C2_Av717Trunc/analysisAlphaV717/furtherAnalysis/adaptiveWindows/figuresSptVsWindows
h(4) = open('fractionMode1_protType1.fig');
cd /home/kj35/files/LCCB/receptors/Galbraiths/data/alphaV717TruncAndCellEdge/120517_Cs2C3_Av717Trunc/analysisAlphaV717/furtherAnalysis/adaptiveWindows/figuresSptVsWindows
h(5) = open('fractionMode1_protType1.fig');
cd /home/kj35/files/LCCB/receptors/Galbraiths/data/alphaV717TruncAndCellEdge/120517_Cs2C1_Av717Trunc/analysisAlphaV717/furtherAnalysis/adaptiveWindows/figuresSptVsWindows
h(6) = open('fractionMode1_protType1.fig');
cd /home/kj35/files/LCCB/receptors/Galbraiths/data/alphaVTruncAndCellEdge/120605_Cs3C3_AvTrunc/analysisAlphaVTrunc/furtherAnalysis/adaptiveWindows/figuresSptVsWindows
h(7) = open('fractionMode1_protType1.fig');
cd /home/kj35/files/LCCB/receptors/Galbraiths/data/alphaVTruncAndCellEdge/120605_Cs1C3_AvTrunc/analysisAlphaVTrunc/furtherAnalysis/adaptiveWindows/figuresSptVsWindows
h(8) = open('fractionMode1_protType1.fig');
cd /home/kj35/files/LCCB/receptors/Galbraiths/data/alphaVTruncAndCellEdge/120217/120217_Cs2C3_AvTrunc/right/analysisAlphaVTrunc/furtherAnalysis/adaptiveWindows/figuresSptVsWindows
h(9) = open('fractionMode1_protType1.fig');
cd /home/kj35/files/LCCB/receptors/Galbraiths/data/alphaVTruncAndCellEdge/120217/120217_Cs2C3_AvTrunc/left/analysisAlphaVTrunc/furtherAnalysis/adaptiveWindows/figuresSptVsWindows
h(10) = open('fractionMode1_protType1.fig');
cd /home/kj35/files/LCCB/receptors/Galbraiths/data/alphaVTruncAndCellEdge/120217/120217_Cs1C3_AvTrunc/analysisAlphaVTrunc/furtherAnalysis/adaptiveWindows/figuresSptVsWindows
h(11) = open('fractionMode1_protType1.fig');
cd /home/kj35/files/LCCB/receptors/Galbraiths/data/alphaVTruncAndCellEdge/120217/120217_Cs1C2_AvTrunc/analysisAlphaVTrunc/furtherAnalysis/adaptiveWindows/figuresSptVsWindows
h(12) = open('fractionMode1_protType1.fig');

for i=1:12
    ha(i) = gca(h(i));
end

for i=1:12
    axis(ha(i),[-3.2 20 0 0.3]);
end


