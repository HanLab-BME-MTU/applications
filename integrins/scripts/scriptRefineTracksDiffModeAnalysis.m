clear all
close all

firstMaskFile = '/home/kj35/files/LCCB/receptors/Galbraiths/data/talinAndCellEdge/120607_Cs2C2_Talin/analysisCellEdgeMod/refined_masks/refined_masks_for_channel_1/refined_mask_mod_120607_Cs2C2_Talin_6minEdgeStack_00001.tif';

cd ../../../analysisTalin/tracks/
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

% [modeParam,expParam] = getDiffModes(tracksFinal,5,0.001,1,10,2,'test');
% modeParam
% 
% save('diffusionModeAnalysis','modeParam','expParam');
% 
% plotPropertySpatialMap2D(tracksFinal,'test',1,5,1,[],[],diffAnalysisRes);
