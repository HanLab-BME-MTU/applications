clear all
close all

firstMaskFile = 'C:\kjData\Galbraiths\data\lifeActAndCellEdge\120224\120224_Cs1C1_LifeAct\120224_Cs1C1a_LifeAct\analysisCellEdgeSmall\SegmentationPackage\refined_masks\refined_masks_for_channel_1\refined_mask_120224Cs1C1a_lifeact_0001.tif';

cd ../../../../analysisLifeact/furtherAnalysis/
load tracksDiffusionLength5InMask.mat

mkdir randomizationTest
cd randomizationTest

seqOfEvents = vertcat(tracksFinal(end-10:end).seqOfEvents);
maxFrame = max(seqOfEvents(:,1));

tracksFinal = randomizeTracksInCellMovie(tracksFinal,firstMaskFile,1:400:maxFrame+1);

save('tracksLength5InMaskRandom.mat','tracksFinal');

