function errArray = frecErrAnalysis(mechDir,simuTrueDir,simuDirList)
%frecErrAnalysis: This function performs the error analysis of force reconstruction on noisy 
%                 simulation data.
%
% SYNOPSIS:
%    errArray = frecErrAnalysis(mechDir,simuTrueDir,simuDirList);
%
% INPUT:
%    mechDir    : The directory of force reconstruction.
%    simuTrueDir: The sub- result directory that contains the true force field used for simulation.
%    simuDirList: A cell array of sub- result directories that contain the reconstructed force
%                 fields from data with various level of noise. Each of such directory corresponds
%                 to one noise level.
%
% OUTPUT:
%    errArray: A structure that contains the error analysis results.
%
%    The result 'errArray' will also be saved under 'mechDir' as 'errArray.mat'.

errArray.mechDir     = mechDir;
errArray.simuTrueDir = simuTrueDir;
errArray.simuDirList = simuDirList;

if ~iscell(simuDirList)
   simuDirList = {simuDirList};
end

%First, load the true force and displacement field.
simuTrueDir = [mechDir filesep simuTrueDir];
forceFieldFile = [simuTrueDir filesep 'forceField' filesep 'forceField.mat'];
s = load(forceFieldFile);
tForceField = s.forceField;

tForceLen          = sqrt(sum(tForceField.f.^2,2));
tinyForceThreshold = nanmean(tForceLen)/2;
iTinyForce         = find(tForceLen<tinyForceThreshold);
iSigForce          = find(tForceLen>=tinyForceThreshold);
meanTrueForceLen   = nanmean(tForceLen(iSigForce));

%True boundary force.
tBndF = tForceField.bndF;
tBndFLen = sqrt(tBndF.fx.^2+tBndF.fy.^2);
tinyBndFThreshold = nanmean(tBndFLen)/10;

answer = input('Edge No. of simulated boundary force:');
if isempty(answer)
   simEdgNo = 1;
else
   simEdgNo = answer;
end

numNoiseLevels = length(simuDirList);

errArray.numNoiseLevels   = numNoiseLevels;
errArray.iTinyForce       = iTinyForce;
errArray.iSigForce        = iSigForce;
errArray.medSpdMean       = zeros(size(numNoiseLevels));
errArray.medSpdStd        = zeros(size(numNoiseLevels));
errArray.absNoiseLevel    = zeros(size(numNoiseLevels));
errArray.relNoiseLevel    = zeros(size(numNoiseLevels));
errArray.bndfRelErrMean   = zeros(size(numNoiseLevels));
errArray.bndfRelErrStd    = zeros(size(numNoiseLevels));
errArray.fRelErrMean      = zeros(size(numNoiseLevels));
errArray.fRelErrStd       = zeros(size(numNoiseLevels));
errArray.fDistrErrMean    = zeros(size(numNoiseLevels));
errArray.fDistrErrStd     = zeros(size(numNoiseLevels));
errArray.errZeroForceMean = zeros(size(numNoiseLevels));
errArray.errZeroForceStd  = zeros(size(numNoiseLevels));
errArray.errZeroFAreaMean = zeros(size(numNoiseLevels));
errArray.errZeroFAreaStd  = zeros(size(numNoiseLevels));
errArray.vRelResMean      = zeros(size(numNoiseLevels));
errArray.vRelResStd       = zeros(size(numNoiseLevels));

for jj = 1:numNoiseLevels
   simuDir = [mechDir filesep simuDirList{jj}];

   [rawDispFieldFiles,index1] = getNamedFiles([simuDir filesep 'rawDispField'],'rawDispField');
   [iDispFieldFiles,index2]   = getNamedFiles([simuDir filesep 'iDispField'],'iDispField');
   [forceFieldFiles,index3]   = getNamedFiles([simuDir filesep 'forceField'],'forceField');
   if length(index1) ~= length(index2) || length(index1) ~= length(index3)
      error(['The number of field files in ''rawDispField'', ' ...
         '''iDispField'' and ''forceField'' do not match.']);
   end
   numRealizations = length(index1);
   s = load([simuDir filesep 'rawDispField' filesep rawDispFieldFiles{1}]);
   errArray.absNoiseLevel(jj) = s.rawDispField.absNoiseLevel;
   errArray.relNoiseLevel(jj) = s.rawDispField.relNoiseLevel;
   medSpd       = [];
   fRelErr      = [];
   fDistrErr    = [];
   bndfRelErr   = [];
   vRelResidue  = [];
   errZeroFArea = [];
   errZeroForce = [];
   for ii = 1:numRealizations
      s = load([simuDir filesep 'rawDispField' filesep rawDispFieldFiles{ii}]);
      rawDispField = s.rawDispField;
      s = load([simuDir filesep 'iDispField' filesep iDispFieldFiles{ii}]);
      iDispField = s.iDispField;
      s = load([simuDir filesep 'forceField' filesep forceFieldFiles{ii}]);
      forceField = s.forceField;

      bndF = forceField.bndF(simEdgNo);

      medSpd  = [medSpd rawDispField.medSpd];
      fRelErr = [fRelErr nanmean(sqrt(sum((forceField.f(iSigForce,:) - ...
                 tForceField.f(iSigForce,:)).^2,2))./tForceLen(iSigForce))];
      fLen    = sqrt(sum(forceField.f.^2,2));
      meanForceLen = nanmean(fLen(iSigForce));
      errZeroFArea = [errZeroFArea ...
                      length(find(fLen(iTinyForce)>tinyForceThreshold))/length(iTinyForce)];
      errZeroForce = [errZeroForce nanmean(sqrt(sum((forceField.f(iTinyForce,:)- ...
                      tForceField.f(iTinyForce,:)).^2,2))./meanTrueForceLen)];
      bndfRelErr   = [bndfRelErr nanmean(sqrt((bndF.fx- ...
                      tBndF.fx).^2 + (bndF.fy-tBndF.fy).^2)./max(tBndFLen,tinyBndFThreshold))];

      dispLen           = sqrt(sum(iDispField.rv.^2,2));
      tinyDispThreshold = nanmean(dispLen)/2;
      vRelResidue       = [vRelResidue nanmean(sqrt(sum((iDispField.v-iDispField.rv).^2,2))./ ...
                           max(tinyDispThreshold,dispLen))];

      %To calculate the error in the reconstructed force distribution, we need to normalize the
      % force field so that the mean of the reconstructed force and the true force is the same. In
      % this way, the error reflects the error in direction, location and relative magnitude while 
      % the difference in the magnitude is ignored.
      fDistrErr = [fDistrErr nanmean(sqrt(sum((forceField.f(iSigForce,:)/meanForceLen - ...
                   tForceField.f(iSigForce,:)/meanTrueForceLen).^2,2))./ ...
                   tForceLen(iSigForce)/meanTrueForceLen)];
      if errArray.absNoiseLevel(jj) ~= rawDispField.absNoiseLevel || ...
         errArray.relNoiseLevel(jj) ~= rawDispField.relNoiseLevel
         error('The noise level should be the same for all realizations.');
      end
   end
   errArray.medSpdMean(jj) = mean(medSpd);
   errArray.medSpdStd(jj)  = std(medSpd);

   %Calculate the relative error in reconstructed forces (>tinyForceThreshold).
   errArray.fRelErrMean(jj)   = mean(fRelErr);
   errArray.fRelErrStd(jj)    = std(fRelErr);
   errArray.fDistrErrMean(jj) = mean(fDistrErr);
   errArray.fDistrErrStd(jj)  = std(fDistrErr);

   %Calculate the area ratio of the reconstructed forces that are greater than 
   % 'tinyForceThreshold' in the tiny force region.
   errArray.errZeroFAreaMean(jj) = mean(errZeroFArea);
   errArray.errZeroFAreaStd(jj)  = std(errZeroFArea);
   errArray.errZeroForceMean(jj) = mean(errZeroForce);
   errArray.errZeroForceStd(jj)  = std(errZeroForce);

   %Calculate the relative error in reconstructed boundary forces (>tinyBndFThreshold).
   errArray.bndfRelErrMean(jj) = mean(bndfRelErr);
   errArray.bndfRelErrStd(jj)  = std(bndfRelErr);

   %Get the residue in matching displacement field.
   errArray.vRelResMean(jj) = mean(vRelResidue);
   errArray.vRelResStd(jj)  = std(vRelResidue);
end

figure; hold on;
plot(1:numNoiseLevels,errArray.bndfRelErrMean,'cs-');
plot(1:numNoiseLevels,errArray.fRelErrMean,'ro-');
plot(1:numNoiseLevels,errArray.errZeroFAreaMean,'gv-');
plot(1:numNoiseLevels,errArray.errZeroForceMean,'k^-');
plot(1:numNoiseLevels,errArray.vRelResMean,'bd-');

errArrayFile = [mechDir filesep 'errArray.mat'];
save(errArrayFile,'errArray');
