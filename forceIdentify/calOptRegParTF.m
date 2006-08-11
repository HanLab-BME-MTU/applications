%For boundary force:
%Identify the optimal regularization parameter range and choose the regularization parameter
% in thisrange.
regParam.minTFSigma = NaN*ones(1,length(imgIndexOfDTimePts));
regParam.maxTFSigma = NaN*ones(1,length(imgIndexOfDTimePts));
regParam.selTFSigma = NaN*ones(1,length(imgIndexOfDTimePts));

regParamFile = [reslDir filesep 'lastSavedRegParam.mat'];
fprintf(1,'Identifying optimal regularization parameter range using L-curve ...\n');
isTFSigmaRangeIdentified = 'no';

answer = input('Select time steps (0 for all):');
if isempty(answer) | answer == 0
   selTimeSteps = 1:numDTimePts;
else
   selTimeSteps = answer;
end

for ii = 1:length(selTimeSteps)
   jj = selTimeSteps(ii);

   fprintf(1,'Time Step %d: ',jj);

   if exist(regParamFile,'file') == 2
      fprintf(1,'\n');
      s = load(regParamFile);
      regParam = s.regParam;

      if ~isfield(regParam,'tfSigma')
         regParam.tfSigma = tfSigma*ones(1,length(imgIndexOfDTimePts));
      end

      if isfield(regParam,'minTFSigma') && isfield(regParam,'maxTFSigma')
         if length(regParam.minTFSigma) == 1
            regParam.minTFSigma = regParam.minTFSigma*ones(1,length(imgIndexOfDTimePts));
            regParam.maxTFSigma = regParam.maxTFSigma*ones(1,length(imgIndexOfDTimePts));
            regParam.selTFSigma = regParam.selTFSigma*ones(1,length(imgIndexOfDTimePts));
         end

         answer = input(sprintf(['   Last saved optimal range of regularization parameter: ' ...
            '(%3.2f,%3.2f)\n' '   Do you want to reidentify it using L-curve? (y/n):'], ...
            regParam.minTFSigma(min(jj,length(regParam.minTFSigma))), ...
            regParam.maxTFSigma(min(jj,length(regParam.minTFSigma)))),'s');
      else
         answer = 'yes';
      end
   end

   %Identify the optimal range of regularization parameter using L-curve.
   if strcmp(answer,'n')
      answer = input(sprintf('   Choose your regularization parameter in this range (last saved: %3.2f):', ...
         regParam.selTFSigma(jj)),'s');
      if strcmp(answer,'')
         selTFSigma = regParam.selTFSigma;
      else
         selTFSigma = str2num(answer);
         regParam.selTFSigma(jj) = selTFSigma;
         save(regParamFile,'regParam');
      end
      continue;
   end

   if isfield(regParam,'tfSigma')
      answer = input(sprintf(['   Please select a new center of test range ' ...
         '(last saved: %5.2f): '],regParam.tfSigma(jj)),'s');
   else
      answer = input('   Please select a new center of test range:','s');
   end

   if ~isempty(answer)
      regParam.tfSigma(jj) = str2num(answer);
   end

   tfSigma = regParam.tfSigma(jj);
   isTFSigmaRangeIdentified = 'no';

   imgIndex = imgIndexOfDTimePts(jj);
   fwdMapTFFile = [fwdMapTFDir filesep 'A' sprintf(imgIndexForm,imgIndex) '.mat'];
   s = load(fwdMapTFFile);
   A = s.A;

   iDispFieldFileName = ['iDispField' sprintf(imgIndexForm,imgIndex) '.mat'];
   iDispFieldFile     = [iDispFieldDir filesep iDispFieldFileName];
   s = load(iDispFieldFile);
   iDispField = s.iDispField;
   rightU     = iDispField.tfRHU;
   numDP      = iDispField.numDP;

   if strcmp(isFieldBndFixed,'yes')
      femModelFile = [femModelDir filesep 'femModel' sprintf(imgIndexForm,0) '.mat'];
   else
      [femModelFile,femModelImgIndex] = getFemModelFile(femModelDir,imgIndex,imgIndexForm);
   end
   s = load(femModelFile);
   femModel = s.femModel;
   fem      = femModel.fem;
   fsBnd    = femModel.fsBnd;

   while strcmp(isTFSigmaRangeIdentified,'no')
      tfSigmaRange2 = tfSigma*[1e-4 1e-3 1e-2 1e-1 1 1e1 1e2 1e3 1e4];
      tfSigmaRange1 = tfSigmaRange2/5;
      tfSigmaRange3 = tfSigmaRange2*5;
      tfSigmaCandidate = [tfSigmaRange1 tfSigmaRange2 tfSigmaRange3];
      tfSigmaCandidate = sort(tfSigmaCandidate);

      coefNorm    = zeros(1,length(tfSigmaCandidate));
      residueNorm = zeros(1,length(tfSigmaCandidate));
      for ik = 1:numEdges
         %Reshape 'A';
         A{ik} = reshape(A{ik},2*numDP,2*fsBnd(ik).dim);
         [m,n] = size(A{ik});

         B0 = A{ik}.'*A{ik};
         b  = A{ik}.'*rightU{ik};

         %Calculate the L-curve.
         for kk = 1:length(tfSigmaCandidate)
            B_i = B0 + tfSigmaCandidate(kk)*eye(n);
            coef_i = B_i\b;
            coefNorm(kk)    = coefNorm(kk)+sqrt(sum(coef_i.^2));
            residueNorm(kk) = residueNorm(kk)+sqrt(sum((A{ik}*coef_i-rightU{ik}).^2));
         end
      end

      %Plot the L-curve.
      figure; hold on;
      plot(coefNorm,residueNorm,'.');
      for kk = 1:length(tfSigmaCandidate)
         text(coefNorm(kk),residueNorm(kk),num2str(tfSigmaCandidate(kk)));
      end

      answer = input('   Do you see a corner of the curve? (yes/no):','s');
      if strcmp(answer,'yes')
         answer = input('   Please identify the corner:','s');
         regParam.selTFSigma(jj) = str2num(answer);
         regParam.minTFSigma(jj) = regParam.selTFSigma(jj)*1e-1;
         regParam.maxTFSigma(jj) = regParam.selTFSigma(jj)*1e1;
         isTFSigmaRangeIdentified = 'yes';
         save(regParamFile,'regParam');
      else
         answer = input('   Please select a new center of test range:','s');
         tfSigma = str2num(answer);
         isTFSigmaRangeIdentified = 'no';
      end
   end
   isTFSigmaRangeIdentified = 'no';

   fprintf(1,'   Optimal regularization parameter range: (%3.2f,%3.2f).\n', ...
      regParam.minTFSigma(jj),regParam.maxTFSigma(jj));
   answer = input(sprintf('Choose your regularization parameter in this range (last saved: %3.2f):', ...
      regParam.selTFSigma(jj)),'s');
   if ~strcmp(answer,'')
      selTFSigma = str2num(answer);
      regParam.selTFSigma(jj) = selTFSigma;
      save(regParamFile,'regParam');
   end
end

isTFSigmaRangeIdentified = 'yes';

