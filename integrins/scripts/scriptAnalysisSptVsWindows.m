clear all
close all

% windowRange = [];
windowRange = (1:135)';
frameRange = [];

% windowsAll = putWindowsTogether;
% load ../../../analysisCellEdgeModSmall/protrusion_samples/protrusion_samples.mat

load ../tracksDiffusionLength5InMask.mat

seqOfEvents = vertcat(tracksFinal(end-10:end).seqOfEvents);
maxFrame = max(seqOfEvents(:,1));

% [windowTrackAssign,trackWindowAssign,trackWindowAssignComp,windowTrackAssignExt] = ...
%     assignTracks2Windows(tracksFinal,windowsAll,1:400:maxFrame+1,1);

% save('windowsActivityTracks','protSamples','windowsAll','trackWindowAssign',...
%     'trackWindowAssignComp','windowTrackAssign','windowTrackAssignExt')

% save('windowsActivityTracks1','protSamples','windowsAll')
% save('windowsActivityTracks2','trackWindowAssign','trackWindowAssignComp','windowTrackAssign')
% size1 = size(windowTrackAssignExt,4);
% size1 = floor(size1/2);
% windowTrackAssignExt1 = windowTrackAssignExt(:,:,:,1:size1);
% windowTrackAssignExt2 = windowTrackAssignExt(:,:,:,size1+1:end);
% save('windowsActivityTracks3','windowTrackAssignExt1')
% save('windowsActivityTracks4','windowTrackAssignExt2')
% clear windowTrackAssignExt1 windowTrackAssignExt2

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
load ../directTrackChar.mat

[sptPropInWindow,~,~,analysisParam] = sptRelToActivityOnsetAdaptiveWindows(...
    tracksFinal,diffAnalysisRes,diffModeAnalysisRes,trackChar,windowsAll,...
    1:400:maxFrame+1,protSamples,windowMSAssignExt,[5 99],windowRange,frameRange,...
    windowTrackAssignExt);

% save('particleBehaviorAdaptiveWindows2','sptPropInWindow','analysisParam');
