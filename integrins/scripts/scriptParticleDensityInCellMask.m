clear all
close all

firstMaskFile = 'C:\kjData\Galbraiths\data\alphaVandCellEdge\120907\120907_Cs2C3\analysisCellEdgeSmall\SegmentationPackage\refined_masks\refined_masks_for_channel_1\refined_mask_120907_Cs2C3_CHO_mEos2Av_6minEdgeStack_0001.tif';

cd ../../../../analysisAlphaV/furtherAnalysis/
load tracksDiffusionLength5InMask.mat

smDensity = particleDensityFromTracksInMask(tracksFinal,firstMaskFile,[5 99]);

save('particleDensity','smDensity');

