clear all
close all

firstMaskFile = '/home/kj35/files/LCCB/receptors/Galbraiths/data/lifeActAndCellEdge/120224/120224_Cs3C1_LifeAct/analysisCellEdgeModSmall/refined_masks/refined_masks_for_channel_1/refined_mask_mod_120224_Cs3C1_CHO_mEos2LifeAct25ms_GFPfill5ms_6minEdgeStack_00001.tif';

cd ../../../analysisLifeact/tracks/
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

