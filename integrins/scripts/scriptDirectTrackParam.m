clear all

load tracksDiffusionLength5InMask.mat
load ../../analysisCellEdgeMod/protrusion_samples/protrusion_samples.mat

windowsAll = putWindowsTogether;

protVec = protSamples.avgVector;
protVecMag = sqrt(sum(protVec.^2,3));
protVecUnit = protVec ./ repmat(protVecMag,[1 1 2]);

seqOfEvents = vertcat(tracksFinal(end-10:end).seqOfEvents);
maxFrame = max(seqOfEvents(:,1));

[windowTrackAssign,trackWindowAssign,trackWindowAssignComp] = assignTracks2Windows(tracksFinal,windowsAll,1:400:maxFrame+1,1);

trackChar = trackMotionCharProtrusion(tracksFinal,protSamples,trackWindowAssignComp,5);

save('directTrackChar','trackChar');
save('windowsActivityTracks','protSamples','windowTrackAssign','trackWindowAssign','trackWindowAssignComp','windowsAll');

