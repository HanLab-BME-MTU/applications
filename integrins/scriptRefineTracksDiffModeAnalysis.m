clear all

maskDir = '/home/kj35/files/LCCB/receptors/Galbraiths/data/alphaVandCellEdge/110114/Cs1_CHO02/Cs1_CHO02A/analysisCellEdgeMod/refined_masks/refined_masks_for_channel_1/refined_mask_mod_110114_Cs1_CHO02A_new_0001.tif';

cd ../../../analysisAlphaV/tracks/
load tracks1AllFrames.mat
cd ../diffusion/
load diffAnalysis1AllFrames.mat
cd ../furtherAnalysis

criteria.lifeTime.min = 5;
indx5 = chooseTracks(tracksFinal,criteria);
tracksFinal = tracksFinal(indx5);
diffAnalysisRes = diffAnalysisRes(indx5);

maxFrame = tracksFinal(end).seqOfEvents(end,1);
indxMask = findTracksInCellMask(tracksFinal,maskDir,1:400:maxFrame+1);
tracksFinal = tracksFinal(indxMask);
diffAnalysisRes = diffAnalysisRes(indxMask);

save('tracksDiffusionLength5InMask','tracksFinal','diffAnalysisRes','indx5','indxMask');

[modeParam,expParam] = getDiffModes(tracksFinal,5,0.001,1,10,2,'test');

save('diffusionModeAnalysis','modeParam','expParam');
