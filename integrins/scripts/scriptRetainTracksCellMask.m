clear all
close all

firstMaskFile = 'C:\kjData\Galbraiths\data\beta3andCellEdgeP2\130517_Cs2C2_Beta3\analysisCellEdgeSmall\SegmentationPackage\refined_masks\refined_masks_for_channel_1\refined_mask_130517_Cs2C2_Beta3_6mES_0001.tif';

cd ../../../../analysisBeta3/tracks/
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
% indxMask = findTracksInCellMask(tracksFinal,firstMaskFile,1:maxFrame:maxFrame+1);
tracksFinal = tracksFinal(indxMask);
diffAnalysisRes = diffAnalysisRes(indxMask);

save('tracksDiffusionLength5InMask','tracksFinal','diffAnalysisRes','indx5','indxMask');

