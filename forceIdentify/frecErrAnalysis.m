%This script file performs the error analysis of force reconstruction on noisy simulation data.

%First, load the true force and displacement field.

forceFieldFile = [forceFieldDir filesep 'forceField.mat'];
s = load(forceFieldFile);
tForceField = s.forceField;
tForceLen          = sqrt(sum(tForceField.f.^2,2));
tinyForceThreshold = nanmean(tForceLen)/2;
iTinyForce         = find(tForceLen<tinyForceThreshold);
iSigForce          = find(tForceLen>=tinyForceThreshold);
meanForceLen       = nanmean(tForceLen(iSigForce));

%True boundary force.
tBndF = tForceField.bndF;
tBndFLen = sqrt(tBndF.fx.^2+tBndF.fy.^2);
tinyBndFThreshold = nanmean(tBndFLen)/10;

answer = input('Select time steps (0 for all):');
if isempty(answer) | answer == 0
   selTimeSteps = 1:numDTimePts;
else
   selTimeSteps = answer;
end

answer = input('Edge No. of simulated boundary force:');
if isempty(answer)
   simEdgNo = 1;
else
   simEdgNo = answer;
end

errArray.iTinyForce    = iTinyForce;
errArray.iSigForce     = iSigForce;
errArray.imgIndex      = imgIndexOfDTimePts(selTimeSteps);
errArray.medSpd        = zeros(size(selTimeSteps));
errArray.absNoiseLevel = zeros(size(selTimeSteps));
errArray.relNoiseLevel = zeros(size(selTimeSteps));
errArray.bndfRelErr    = zeros(size(selTimeSteps));
errArray.fRelErr       = zeros(size(selTimeSteps));
errArray.errZeroForce  = zeros(size(selTimeSteps));
errArray.errZeroFArea  = zeros(size(selTimeSteps));
errArray.vRelRes       = zeros(size(selTimeSteps));
for ii = 1:length(selTimeSteps)
   jj = selTimeSteps(ii);

   imgIndex = imgIndexOfDTimePts(jj);

   rawDispFieldFileName = ['rawDispField' sprintf(imgIndexForm,imgIndex) '.mat'];
   rawDispFieldFile     = [rawDispFieldDir filesep rawDispFieldFileName];
   s = load(rawDispFieldFile);
   rawDispField = s.rawDispField;
   errArray.medSpd(ii)        = rawDispField.medSpd;
   errArray.absNoiseLevel(ii) = rawDispField.absNoiseLevel;
   errArray.relNoiseLevel(ii) = rawDispField.relNoiseLevel;

   iDispFieldFileName = ['iDispField' sprintf(imgIndexForm,imgIndex) '.mat'];
   iDispFieldFile     = [iDispFieldDir filesep iDispFieldFileName];
   s = load(iDispFieldFile);
   iDispField = s.iDispField;

   forceFieldFile = [forceFieldDir filesep 'forceField' ...
      sprintf(imgIndexForm,imgIndex) '.mat'];
   s = load(forceFieldFile);
   forceField = s.forceField;
   bndF       = forceField.bndF(simEdgNo);

   %Calculate the relative error in reconstructed forces (>tinyForceThreshold).
   errArray.fRelErr(ii) = nanmean(abs(sqrt(sum((forceField.f(iSigForce,:) - ...
      tForceField.f(iSigForce,:)).^2,2))./tForceLen(iSigForce)));

   %Calculate the area ratio of the reconstructed forces that are greater than 
   % 'tinyForceThreshold' in the tiny force region.
   fLen = sqrt(sum(forceField.f.^2,2));
   errArray.errZeroFArea(ii) = length(find(fLen(iTinyForce)>tinyForceThreshold))/length(iTinyForce);
   errArray.errZeroForce(ii) = nanmean(sqrt(sum((forceField.f(iTinyForce,:)- ...
      tForceField.f(iTinyForce,:)).^2,2))./meanForceLen);

   %Calculate the relative error in reconstructed boundary forces (>tinyBndFThreshold).
   errArray.bndfRelErr(ii) = nanmean(sqrt((bndF.fx- ...
      tBndF.fx).^2 + (bndF.fy-tBndF.fy).^2)./max(tBndFLen,tinyBndFThreshold));

   %Get the residue in matching displacement field.
   dispLen          = sqrt(sum(iDispField.rv.^2,2));
   tinyDispThreshold = nanmean(dispLen)/2;
   errArray.vRelRes(ii) = nanmean(abs(sqrt(sum((iDispField.v-iDispField.rv).^2,2))./ ...
      max(tinyDispThreshold,dispLen)));
end

figure; hold on;
plot(errArray.imgIndex,errArray.bndfRelErr,'cs-');
plot(errArray.imgIndex,errArray.fRelErr,'ro-');
plot(errArray.imgIndex,errArray.errZeroFArea,'gv-');
plot(errArray.imgIndex,errArray.errZeroForce,'k^-');
plot(errArray.imgIndex,errArray.vRelRes,'bd-');

errArrayFile = [reslDir filesep 'errArray.mat'];
save(errArrayFile,'errArray');
