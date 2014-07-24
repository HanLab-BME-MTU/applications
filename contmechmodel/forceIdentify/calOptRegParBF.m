%Identify the optimal regularization parameter range and choose the regularization parameter
% in thisrange.

%regParam.minBFSigma = NaN*ones(1,length(imgIndexOfDTimePts));
%regParam.maxBFSigma = NaN*ones(1,length(imgIndexOfDTimePts));
%regParam.selBFSigma = NaN*ones(1,length(imgIndexOfDTimePts));

clear regParam;

regParamFile = [reslDir filesep 'lastSavedRegParam.mat'];
fprintf(1,'Identifying optimal regularization parameter range using L-curve ...\n');
isBFSigmaRangeIdentified = 'no';

answer = input('Select time steps (0 for all):');
if isempty(answer) | answer == 0
   selTimeSteps = 1:numDTimePts;
else
   selTimeSteps = answer;
end

if exist(regParamFile,'file') == 2
   s = load(regParamFile);
   regParam = s.regParam;
end

for ii = 1:length(selTimeSteps)
   jj = selTimeSteps(ii);

   fprintf(1,'Time Step %d: ',jj);
   imgIndex = imgIndexOfDTimePts(jj);

   clear forceField;

   forceFieldFile = [forceFieldDir filesep 'forceField' ...
      sprintf(imgIndexForm,imgIndex) '.mat'];

   selBFSigma = [];
   if exist(forceFieldFile,'file') == 2
      s = load(forceFieldFile);
      forceField = s.forceField;

      if isfield(forceField,'selBFSigma')
         bfSigma    = forceField.bfSigma;
         minBFSigma = forceField.minBFSigma;
         maxBFSigma = forceField.maxBFSigma;
         selBFSigma = forceField.selBFSigma;

         %Synchronize between 'forceField' and 'regParam'.
         regParam.bfSigma(jj)    = bfSigma;
         regParam.selBFSigma(jj) = selBFSigma;
         regParam.minBFSigma(jj) = minBFSigma;
         regParam.maxBFSigma(jj) = maxBFSigma;

         save(regParamFile,'regParam');
      end
   end

   if isempty(selBFSigma)
      if exist(regParamFile,'file') == 2
         fprintf(1,'\n');
         s = load(regParamFile);
         regParam = s.regParam;

         if isfield(regParam,'bfSigma')
            bfSigma = regParam.bfSigma(min(jj,length(regParam.bfSigma)));
         end
         minBFSigma = regParam.minBFSigma(min(jj,length(regParam.minBFSigma)));
         maxBFSigma = regParam.maxBFSigma(min(jj,length(regParam.maxBFSigma)));
         selBFSigma = regParam.selBFSigma(min(jj,length(regParam.selBFSigma)));

         %Synchronize between 'forceField' and 'regParam'.
         regParam.bfSigma(jj)    = bfSigma;
         regParam.selBFSigma(jj) = selBFSigma;
         regParam.minBFSigma(jj) = minBFSigma;
         regParam.maxBFSigma(jj) = maxBFSigma;

         forceField.bfSigma    = bfSigma;
         forceField.minBFSigma = minBFSigma;
         forceField.maxBFSigma = maxBFSigma;
         forceField.selBFSigma = selBFSigma;

         save(regParamFile,'regParam');
         save(forceFieldFile,'forceField');
      end
   end

   if ~isempty(selBFSigma)
      answer = input(sprintf(['   Last saved optimal range of regularization parameter: ' ...
         '(%3.2f,%3.2f)\n' '   Do you want to reidentify it using L-curve? (y/n):'], ...
         minBFSigma,maxBFSigma),'s');
   else
      answer = 'y';
   end

   %Identify the optimal range of regularization parameter using L-curve.
   if strcmp(answer,'n')
      regParam.bfSigma(jj)    = bfSigma;
      regParam.minBFSigma(jj) = minBFSigma;
      regParam.maxBFSigma(jj) = maxBFSigma;

      forceField.bfSigma    = bfSigma;
      forceField.minBFSigma = minBFSigma;
      forceField.maxBFSigma = maxBFSigma;

      answer = input(sprintf('   Choose your regularization parameter in this range (last saved: %3.2f):', ...
         selBFSigma),'s');
      if ~strcmp(answer,'')
         selBFSigma = str2num(answer);
      end

      regParam.selBFSigma(jj) = selBFSigma;
      forceField.selBFSigma   = selBFSigma;

      save(regParamFile,'regParam');
      save(forceFieldFile,'forceField');
      continue;
   end

   answer = input(sprintf(['   Please select a new center of test range ' ...
      '(last saved: %5.2f): '],bfSigma),'s');

   if ~isempty(answer)
      bfSigma = str2num(answer);
   end

   regParam.bfSigma(jj) = bfSigma;
   forceField.bfSigma   = bfSigma;

   save(regParamFile,'regParam');
   save(forceFieldFile,'forceField');
   isBFSigmaRangeIdentified = 'no';

   if strcmp(isFieldBndFixed,'yes') && strcmp(isDataPosFixed,'yes')
      fwdMapBFFile = [fwdMapBFDir filesep 'A' sprintf(imgIndexForm,0) '.mat'];
   else
      fwdMapBFFile = [fwdMapBFDir filesep 'A' sprintf(imgIndexForm,imgIndex) '.mat'];
   end
   s = load(fwdMapBFFile);
   A = s.A;

   iDispFieldFileName = ['iDispField' sprintf(imgIndexForm,imgIndex) '.mat'];
   iDispFieldFile     = [iDispFieldDir filesep iDispFieldFileName];
   s = load(iDispFieldFile);
   iDispField = s.iDispField;

   if strcmp(isFieldBndFixed,'yes')
      femModelFile = [femModelDir filesep 'femModel' sprintf(imgIndexForm,0) '.mat'];
   else
      [femModelFile,femModelImgIndex] = getFemModelFile(femModelDir,imgIndex,imgIndexForm);
   end
   s = load(femModelFile);
   femModel = s.femModel;
   fem      = femModel.fem;
   fs       = femModel.fs;

   dimBF     = fs.dimBF;
   dimFS     = fs.dimFS;
   indDomDOF = fs.indDomDOF;

   rightU = iDispField.bfRHU;

   [m,n] = size(A);

   B0 = A.'*A;
   b  = A.'*rightU;

   while strcmp(isBFSigmaRangeIdentified,'no')
      %Calculate the L-curve.
      bfSigmaRange2 = bfSigma*[1e-5 1e-4 1e-3 1e-2 1e-1 1 1e1 1e2 1e3 1e4 1e5];
      bfSigmaRange1 = bfSigmaRange2/5;
      bfSigmaRange3 = bfSigmaRange2*5;
      bfSigmaCandidate = [bfSigmaRange1 bfSigmaRange2 bfSigmaRange3];
      bfSigmaCandidate = sort(bfSigmaCandidate);
      coefNorm    = zeros(1,length(bfSigmaCandidate));
      residueNorm = zeros(1,length(bfSigmaCandidate));
      for kk = 1:length(bfSigmaCandidate)
         B_i = B0 + bfSigmaCandidate(kk)*eye(n);
         coef_i = B_i\b;
         coefNorm(kk)    = sqrt(sum(coef_i.^2));
         residueNorm(kk) = sqrt(sum((A*coef_i-rightU).^2));
      end

      %Plot the L-curve.
      figure; hold on;
      plot(coefNorm,residueNorm,'.');
      for kk = 1:length(bfSigmaCandidate)
         text(coefNorm(kk),residueNorm(kk),num2str(bfSigmaCandidate(kk)));
      end

      answer = input('   Do you see a corner of the curve? (yes/no):','s');
      if strcmp(answer,'yes')
         answer = input('   Please identify the corner:','s');
         regParam.bfSigma(jj)    = bfSigma;
         regParam.selBFSigma(jj) = str2num(answer);
         regParam.minBFSigma(jj) = regParam.selBFSigma(jj)*1e-1;
         regParam.maxBFSigma(jj) = regParam.selBFSigma(jj)*1e1;

         forceField.bfSigma    = bfSigma;
         forceField.selBFSigma = str2num(answer);
         forceField.minBFSigma = regParam.selBFSigma(jj)*1e-1;
         forceField.maxBFSigma = regParam.selBFSigma(jj)*1e1;

         isBFSigmaRangeIdentified = 'yes';
         save(regParamFile,'regParam');
         save(forceFieldFile,'forceField');
      else
         answer = input(sprintf(['   Please select a new center of test range ' ...
            '(last used: %5.2f): '],bfSigma),'s');
         bfSigma = str2num(answer);
         isBFSigmaRangeIdentified = 'no';
      end
   end

   isBFSigmaRangeIdentified = 'no';

   fprintf(1,'   Optimal regularization parameter range: (%3.2f,%3.2f).\n', ...
      forceField.minBFSigma,forceField.maxBFSigma);
   answer = input(sprintf('   Choose your regularization parameter in this range (last saved: %3.2f):', ...
      forceField.selBFSigma),'s');
   if ~strcmp(answer,'')
      selBFSigma = str2num(answer);
      regParam.selBFSigma(jj) = selBFSigma;
      forceField.selBFSigma   = selBFSigma;
      save(regParamFile,'regParam');
      save(forceFieldFile,'forceField');
   end
end

isBFSigmaRangeIdentified = 'yes';

