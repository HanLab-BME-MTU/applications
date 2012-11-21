clear all
close all

firstMaskFile = '/home/kj35/files/LCCB/receptors/Galbraiths/data/alphaV717TruncAndCellEdge/121011_Cs3C2_Av717Trunc/analysisCellEdgeSmall/SegmentationPackage/refined_masks/refined_masks_for_channel_1/refined_mask_121011_Cs3C2_CHO_mEos2Av717Trunc_6minEdgeStack_0001.tif';

cd ../../../../analysisAlphaV717Trunc/tracks/
load tracks1AllFrames.mat
cd ../diffusion/
load diffAnalysis1AllFrames.mat
cd ../furtherAnalysis

criteria.lifeTime.min = 5;
indx5 = chooseTracks(tracksFinal,criteria);
tracksFinal = tracksFinal(indx5);
diffAnalysisRes = diffAnalysisRes(indx5);

seqOfEvents = vertcat(tracksFinal(end-10:end).seqOfEvents);
maxFrame = max(seqOfEvents(:,1));
indxMask = findTracksInCellMask(tracksFinal,firstMaskFile,1:400:maxFrame+1);
tracksFinal = tracksFinal(indxMask);
diffAnalysisRes = diffAnalysisRes(indxMask);

save('tracksDiffusionLength5InMask','tracksFinal','diffAnalysisRes','indx5','indxMask');

