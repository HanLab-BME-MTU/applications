clear all
close all

% windowRange = [];
windowRange = (1:125)';
frameRange = [];

windowsAll = putWindowsTogether;
load ../../../analysisCellEdgeModSmall/protrusion_samples/protrusion_samples.mat

load ../tracksDiffusionLength5InMask.mat

seqOfEvents = vertcat(tracksFinal(end-10:end).seqOfEvents);
maxFrame = max(seqOfEvents(:,1));

[windowTrackAssign,trackWindowAssign,trackWindowAssignComp,windowTrackAssignExt] = ...
    assignTracks2Windows(tracksFinal,windowsAll,1:400:maxFrame+1,1);
save('windowsActivityTracks','protSamples','windowsAll','trackWindowAssign',...
    'trackWindowAssignComp','windowTrackAssign','windowTrackAssignExt')

% load windowsActivityTracks.mat

load ../diffusionModeClassification.mat
load ../directTrackChar.mat

[sptPropInWindow,~,~,analysisParam] = sptRelToActivityOnsetAdaptiveWindows(...
    tracksFinal,diffAnalysisRes,diffModeAnalysisRes,trackChar,windowsAll,...
    1:400:maxFrame+1,protSamples,5,windowRange,frameRange,windowTrackAssignExt);

save('particleBehaviorAdaptiveWindows','sptPropInWindow','analysisParam');
