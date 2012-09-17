


% sliceRange = [];
sliceRange = (40:170)';
frameRange = [];
lengthMinMax = [5 99];

load ../tracksDiffusionLength5InMask.mat

try
    load windowsActivityTracks.mat
catch
    load windowsActivityTracks1.mat
    load windowsActivityTracks2.mat
    load windowsActivityTracks3.mat
    load windowsActivityTracks4.mat
    windowTrackAssignExt = cat(4,windowTrackAssignExt1,windowTrackAssignExt2);
    clear windowTrackAssignExt1 windowTrackAssignExt2
end

load ../diffusionModeClassification.mat
load directTrackChar.mat

load windowNumbersAssignExt.mat

firstMaskFile = '/home/kj35/files/LCCB/receptors/Galbraiths/data/alphaVandCellEdge/110114/Cs1_CHO02/Cs1_CHO02A/analysisCellEdgeSmall2/refined_masks/refined_masks_for_channel_1/refined_mask_110114_Cs1_CHO02A_new_0001.tif';

seqOfEvents = vertcat(tracksFinal(end-10:end).seqOfEvents);
maxFrame = max(seqOfEvents(:,1));

protWinParam = struct('numTypeProt',9,'numPixInBand',2,...
    'numSmallBandsInBigBand',3,'numBigBands',10,'maxNegInc',3,'maxPosInc',18);

[sptPropInWindow,windowDistFromEdge,analysisParam] = sptRelToActivityOnsetAdaptiveWindows(...
    tracksFinal,diffAnalysisRes,diffModeAnalysisRes,trackChar,windowsAll,...
    protSamples,windowTrackAssignExt,windowNumbersAssignExt,...
    lengthMinMax,sliceRange,frameRange,firstMaskFile,protWinParam);

save('particleBehaviorAdaptiveWindows120917','sptPropInWindow',...
    'windowDistFromEdge','analysisParam');

