clear all
close all

mkdir adaptiveWindows
cd adaptiveWindows

% sliceRange = [];
sliceRange = (40:170)';
frameRange = [];
lengthMinMax = [5 99];

windowsAll = putWindowsTogether;
% load ../../../analysisCellEdgeModSmall/protrusion_samples/protrusion_samples.mat
load ../../../analysisCellEdgeSmall2/protrusion_samples/protrusion_samples.mat

load ../tracksDiffusionLength5InMask.mat

seqOfEvents = vertcat(tracksFinal(end-10:end).seqOfEvents);
maxFrame = max(seqOfEvents(:,1));

[windowTrackAssign,trackWindowAssign,trackWindowAssignComp,windowTrackAssignExt] = ...
    assignTracks2Windows(tracksFinal,windowsAll,1:400:maxFrame+1,1);

save('windowsActivityTracks','protSamples','windowsAll','trackWindowAssign',...
    'trackWindowAssignComp','windowTrackAssign','windowTrackAssignExt')

% save('windowsActivityTracks1','protSamples','windowsAll')
% save('windowsActivityTracks2','trackWindowAssign','trackWindowAssignComp','windowTrackAssign')
% size1 = size(windowTrackAssignExt,4);
% size1 = floor(size1/2);
% windowTrackAssignExt1 = windowTrackAssignExt(:,:,:,1:size1);
% windowTrackAssignExt2 = windowTrackAssignExt(:,:,:,size1+1:end);
% save('windowsActivityTracks3','windowTrackAssignExt1')
% save('windowsActivityTracks4','windowTrackAssignExt2')
% clear windowTrackAssignExt1 windowTrackAssignExt2
