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

tForceLen           = sqrt(sum(tForceField.f.^2,2));
tTinyForceThreshold = nanmean(tForceLen)/2;
iTrueTinyForce      = find(tForceLen<tTinyForceThreshold);
iTrueSigForce       = find(tForceLen>=tTinyForceThreshold);
meanTrueForceLen    = nanmean(tForceLen(iTrueSigForce));
tinyForcePercent    = tTinyForceThreshold/meanTrueForceLen;

%True boundary force.
tBndF = tForceField.bndF;
tBndFLen = sqrt(tBndF.fx.^2+tBndF.fy.^2);
meanTrueBndFLen    = nanmean(tBndFLen);
tTinyBndFThreshold = meanTrueBndFLen/10;
iTrueSigBndF       = find(tBndFLen>=tTinyBndFThreshold);

answer = input('Edge No. of simulated boundary force:');
if isempty(answer)
   simEdgNo = 1;
else
   simEdgNo = answer;
end

numNoiseLevels = length(simuDirList);

errArray.numNoiseLevels      = numNoiseLevels;
errArray.tTinyForceThreshold = tTinyForceThreshold;
errArray.iTrueTinyForce      = iTrueTinyForce;
errArray.tinyForcePercent    = tinyForcePercent;
errArray.iTrueSigForce       = iTrueSigForce;
errArray.meanTrueForceLen    = meanTrueForceLen;
errArray.tTinyBndFThreshold  = tTinyBndFThreshold;
errArray.meanTrueBndFLen     = meanTrueBndFLen;
errArray.tinyForceThrMean    = zeros(1,numNoiseLevels);
errArray.tinyForceThrStd     = zeros(1,numNoiseLevels);
errArray.iSigTinyF           = cell(5,numNoiseLevels);
errArray.iTinySigF           = cell(5,numNoiseLevels);
errArray.meanForceLenMean    = zeros(1,numNoiseLevels);
errArray.meanForceLenStd     = zeros(1,numNoiseLevels);
errArray.meanBndFLenMean     = zeros(1,numNoiseLevels);
errArray.meanBndFLenStd      = zeros(1,numNoiseLevels);
errArray.medSpdMean          = zeros(1,numNoiseLevels);
errArray.medSpdStd           = zeros(1,numNoiseLevels);
errArray.absNoiseLevel       = zeros(1,numNoiseLevels);
errArray.relNoiseLevel       = zeros(1,numNoiseLevels);
errArray.bndfRelErrMean      = zeros(1,numNoiseLevels);
errArray.bndfRelErrStd       = zeros(1,numNoiseLevels);
errArray.bndfRelMagErrMean   = zeros(1,numNoiseLevels);
errArray.bndfRelMagErrStd    = zeros(1,numNoiseLevels);
errArray.fRelErrMean         = zeros(1,numNoiseLevels);
errArray.fRelErrStd          = zeros(1,numNoiseLevels);
errArray.fRelMagErrMean      = zeros(1,numNoiseLevels);
errArray.iFRelErrMean        = zeros(1,numNoiseLevels);
errArray.fRelMagErrStd       = zeros(1,numNoiseLevels);
errArray.errZeroForceMean    = zeros(1,numNoiseLevels);
errArray.errZeroForceStd     = zeros(1,numNoiseLevels);
errArray.errZeroFAreaMean    = zeros(5,numNoiseLevels);
errArray.errZeroFAreaStd     = zeros(5,numNoiseLevels);
errArray.vRelResMean         = zeros(1,numNoiseLevels);
errArray.vRelResStd          = zeros(1,numNoiseLevels);

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
   medSpd             = zeros(1,numRealizations);
   iSigTinyF          = cell(5,numRealizations);
   iTinySigF          = cell(5,numRealizations);
   meanForceLen       = zeros(1,numRealizations);
   meanBndFLen        = zeros(1,numRealizations);
   tinyForceThreshold = zeros(1,numRealizations);
   fRelErr            = zeros(1,numRealizations);
   fRelMagErr         = zeros(1,numRealizations);
   fDirErr            = zeros(1,numRealizations);
   bndfRelErr         = zeros(1,numRealizations);
   bndfRelMagErr      = zeros(1,numRealizations);
   bndfDirErr         = zeros(1,numRealizations);
   vRelResidue        = zeros(1,numRealizations);
   errZeroFArea       = zeros(5,numRealizations);
   errZeroForce       = zeros(1,numRealizations);
   for ii = 1:numRealizations
      s = load([simuDir filesep 'rawDispField' filesep rawDispFieldFiles{ii}]);
      rawDispField = s.rawDispField;
      s = load([simuDir filesep 'iDispField' filesep iDispFieldFiles{ii}]);
      iDispField = s.iDispField;
      s = load([simuDir filesep 'forceField' filesep forceFieldFiles{ii}]);
      forceField = s.forceField;

      bndF            = forceField.bndF(simEdgNo);
      bndFLen         = sqrt(bndF.fx.^2+bndF.fy.^2);
      meanBndFLen(ii) = nanmean(bndFLen);

      medSpd(ii)  = rawDispField.medSpd;
      fRelErr(ii) = nanmean(sqrt(sum((forceField.f(iTrueSigForce,:) - ...
         tForceField.f(iTrueSigForce,:)).^2,2))./tForceLen(iTrueSigForce));
      fLen = sqrt(sum(forceField.f.^2,2));

      meanForceLen(ii)       = nanmean(fLen(iTrueSigForce));
      tinyForceThreshold(ii) = meanForceLen(ii)*errArray.tinyForcePercent;

      %We calculate the error in identified zero force area based on 5 levels of tiny force
      % threshold.
      for kk = 1:5
         %Tiny force (below threshold) in the supposed zero force area.
         iSigTinyF{kk,ii} = iTrueTinyForce(find(fLen(iTrueTinyForce)>kk*tinyForceThreshold(ii)));

         %Tiny force (below threshold) in the supposed significant force area.
         iTinySigF{kk,ii} = iTrueSigForce(find(fLen(iTrueSigForce)<=kk*tinyForceThreshold(ii)));

         %The error of identified force location is a combination of the misidentified 
         % significant force area in the supposed zero force region and the misidentified 
         % tiny force area in the supposed signifcant force region.
         errZeroFArea(kk,ii) = (length(iTinySigF{kk,ii})+length(iSigTinyF{kk,ii}))/ ...
            (length(iTrueSigForce)+length(iTrueTinyForce));
      end
      errZeroForce(ii) = nanmean(sqrt(sum((forceField.f(iTrueTinyForce,:)- ...
         tForceField.f(iTrueTinyForce,:)).^2,2))./meanTrueForceLen);
      bndfRelErr(ii)   = nanmean(sqrt((bndF.fx- ...
         tBndF.fx).^2 + (bndF.fy-tBndF.fy).^2)./max(tBndFLen,tTinyBndFThreshold));

      dispLen           = sqrt(sum(iDispField.rv.^2,2));
      tinyDispThreshold = nanmean(dispLen)/2;
      vRelResidue(ii)   = nanmean(sqrt(sum((iDispField.v-iDispField.rv).^2,2))./ ...
         max(tinyDispThreshold,dispLen));

      %To calculate the error in the reconstructed force distribution, we need to normalize the
      % force field so that the mean of the reconstructed force and the true force is the same. In
      % this way, the error reflects the error in direction and relative magnitude while 
      % the difference in the absolute magnitude is ignored.
      %Magnitude:
      fRelMagErr(ii) = nanmean((fLen(iTrueSigForce)/meanForceLen(ii) - ...
         tForceLen(iTrueSigForce)/meanTrueForceLen)./ ...
         (tForceLen(iTrueSigForce)/meanTrueForceLen));
      %Direction:
      fDirErr(ii) = nanmean(abs(acos(sum(forceField.f(iTrueSigForce,:).* ...
         tForceField.f(iTrueSigForce,:),2)./fLen(iTrueSigForce)./tForceLen(iTrueSigForce)))/pi);

      %Relative error in boundary force relative magnitude and direction.
      %Magnitude:
      bndfRelMagErr(ii) = nanmean((bndFLen/meanBndFLen(ii) - ...
         tBndFLen/meanTrueBndFLen)./ ...
         (max(tBndFLen,tTinyBndFThreshold)/meanTrueBndFLen));
      %Direction:
      bndfDirErr(ii) = nanmean(abs(acos((bndF.fx(iTrueSigBndF).*tBndF.fx(iTrueSigBndF) + ...
         bndF.fy(iTrueSigBndF).*tBndF.fy(iTrueSigBndF))./bndFLen(iTrueSigBndF)./ ...
         tBndFLen(iTrueSigBndF)))/pi);
      if errArray.absNoiseLevel(jj) ~= rawDispField.absNoiseLevel || ...
         errArray.relNoiseLevel(jj) ~= rawDispField.relNoiseLevel
         error('The noise level should be the same for all realizations.');
      end
   end
   errArray.medSpdMean(jj) = mean(medSpd);
   errArray.medSpdStd(jj)  = std(medSpd);

   %Calculate the relative error in reconstructed forces in nonzero force region i.e. area where the
   %true forces are > tTinyForceThreshold.
   errArray.fRelErrMean(jj)   = mean(fRelErr);
   errArray.fRelErrStd(jj)    = std(fRelErr);
   errArray.tinyForceThrMean(jj) = mean(tinyForceThreshold);
   errArray.tinyForceThrStd(jj)  = std(tinyForceThreshold);
   errArray.meanForceLenMean(jj) = mean(meanForceLen);
   errArray.meanForceLenStd(jj)  = std(meanForceLen);
   errArray.fRelMagErrMean(jj)   = mean(fRelMagErr);
   errArray.fRelMagErrStd(jj)    = std(fRelMagErr);
   errArray.fDirErrMean(jj)      = mean(fDirErr);
   errArray.fDirErrStd(jj)       = std(fDirErr);

   %Calculate misidentified force area that is defined by area of reconstructed forces that are 
   % greater than 'tinyForceThreshold' in the tiny force region.
   errArray.errZeroFAreaMean(:,jj) = mean(errZeroFArea,2);
   errArray.errZeroFAreaStd(:,jj)  = std(errZeroFArea,0,2);
   errArray.errZeroForceMean(jj)   = mean(errZeroForce);
   errArray.errZeroForceStd(jj)    = std(errZeroForce);

   %Calculate the relative error in reconstructed boundary forces.
   errArray.bndfRelErrMean(jj) = mean(bndfRelErr);
   errArray.bndfRelErrStd(jj)  = std(bndfRelErr);

   %Calculate the relative error in reconstructed boundary force distribution.
   errArray.meanBndFLenMean(jj)   = mean(meanBndFLen);
   errArray.meanBndFLenStd(jj)    = std(meanBndFLen);
   errArray.bndfRelMagErrMean(jj) = mean(bndfRelMagErr);
   errArray.bndfRelMagErrStd(jj)  = std(bndfRelMagErr);
   errArray.bndfDirErrMean(jj)    = mean(bndfDirErr);
   errArray.bndfDirErrStd(jj)     = std(bndfDirErr);

   %Get the residue in matching displacement field.
   errArray.vRelResMean(jj) = mean(vRelResidue);
   errArray.vRelResStd(jj)  = std(vRelResidue);

   %Get the realization whose corresponding error in reconstructed force distribution is closest to
   %the mean error.
   [tmpVar errArray.iFRelErrMean(jj)] = min(abs(fRelErr-errArray.fRelErrMean(jj)));
   errArray.iSigTinyF(:,jj) = iSigTinyF(:,errArray.iFRelErrMean(jj));
   errArray.iTinySigF(:,jj)  = iTinySigF(:,errArray.iFRelErrMean(jj));
end

vRelRes = errArray.vRelResMean;

figure; hold on;
%L11 = errorbar(1:numNoiseLevels,errArray.fRelErrMean,errArray.fRelErrStd);
%L12 = errorbar(1:numNoiseLevels,errArray.bndfRelErrMean,errArray.bndfRelErrStd);
%L13 = errorbar(1:numNoiseLevels,errArray.errZeroForceMean,errArray.errZeroForceStd);
L11 = errorbar(vRelRes,errArray.fRelErrMean,errArray.fRelErrStd);
L12 = errorbar(vRelRes,errArray.bndfRelErrMean,errArray.bndfRelErrStd);
L13 = errorbar(vRelRes,errArray.errZeroForceMean,errArray.errZeroForceStd);
set(L11,'Color','r','LineWidth',5);
set(L11(2),'LineWidth',2);
set(L12,'Color','c','LineWidth',3);
set(L12(2),'LineWidth',2);
set(L13,'Color','b','LineWidth',1);
set(L13(2),'LineWidth',2);
title('Error of Whole Force Vector');
xlabel('Flow Vector Residue');
ylabel('Relative Error');
set(gca,'XTick',vRelRes);

figure; hold on;
L21 = errorbar(vRelRes,errArray.fRelMagErrMean,errArray.fRelMagErrStd);
L22 = errorbar(vRelRes,errArray.bndfRelMagErrMean,errArray.bndfRelMagErrStd);
L23 = errorbar(vRelRes,errArray.fDirErrMean,errArray.fDirErrStd);
L24 = errorbar(vRelRes,errArray.bndfDirErrMean,errArray.bndfDirErrStd);
%L22 = plot(1:numNoiseLevels,errArray.vRelResMean);
set(L21,'Color','r','LineWidth',3);
set(L21(2),'LineWidth',2);
set(L22,'Color','c','LineWidth',3);
set(L22(2),'LineWidth',2);
set(L23,'Color','r','LineWidth',1);
set(L23(2),'LineWidth',2);
set(L24,'Color','c','LineWidth',1);
set(L24(2),'LineWidth',2);
title('Error of Relative Magnitude and Dirction');
xlabel('Flow Vector Residue');
ylabel('Relative Error');
set(gca,'XTick',vRelRes);

figure; hold on;
L3 = errorbar(vRelRes.'*ones(1,5),errArray.errZeroFAreaMean.',errArray.errZeroFAreaStd.');
for kk = 1:5
   set(L3(5-kk+1),'LineWidth',kk);
end
set(L3(1),'Color',[0.8,0.8,0.8]);
set(L3(2),'Color',[0.7,0.7,0.7]);
set(L3(3),'Color',[0.5,0.5,0.5]);
set(L3(4),'Color',[0,0,0]);
set(L3(5),'Color','b','LineWidth',1);
set(L3(6:end),'Color',[0.6,0.6,0.6],'LineWidth',2);
set(L3(end),'Color','b','LineWidth',2);
title('Error of Force Location');
xlabel('Flow Vector Residue');
ylabel('Relative Error');
set(gca,'XTick',vRelRes);
%set(L22,'Color','b','LineWidth',3);

errArrayFile = [mechDir filesep 'errArray.mat'];
save(errArrayFile,'errArray');
