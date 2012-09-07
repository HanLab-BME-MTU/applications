


sliceRange = [];
% sliceRange = (141:201)';
frameRange = [];
lengthMinMax = [5 99];

% load ../tracksDiffusionLength5InMask.mat
% 
% try
%     load windowsActivityTracks.mat
% catch
%     load windowsActivityTracks1.mat
%     load windowsActivityTracks2.mat
%     load windowsActivityTracks3.mat
%     load windowsActivityTracks4.mat
%     windowTrackAssignExt = cat(4,windowTrackAssignExt1,windowTrackAssignExt2);
%     clear windowTrackAssignExt1 windowTrackAssignExt2
% end
% 
% load ../diffusionModeClassification.mat
% load directTrackChar.mat
% 
% load windowNumbersAssignExt.mat

firstMaskFile = '/home/kj35/files/LCCB/receptors/Galbraiths/data/alphaVandCellEdge/110114/Cs1_CHO03/Cs1_CHO03A/analysisCellEdgeModSmall/refined_masks/refined_masks_for_channel_1/refined_mask_mod_110114_Cs1_CHO03A_new_0001.tif';

seqOfEvents = vertcat(tracksFinal(end-10:end).seqOfEvents);
maxFrame = max(seqOfEvents(:,1));

[sptPropInWindow,windowDistFromEdge,analysisParam] = sptRelToActivityOnsetAdaptiveWindows(...
    tracksFinal,diffAnalysisRes,diffModeAnalysisRes,trackChar,windowsAll,...
    protSamples,windowTrackAssignExt,windowNumbersAssignExt,...
    lengthMinMax,sliceRange,frameRange,firstMaskFile);

Sn save('particleBehaviorAdaptiveWindows120907','sptPropInWindow',...
    'windowDistFromEdge','analysisParam');

