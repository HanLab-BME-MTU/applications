%This script file solves the linear system assembled by 'calFwdOpBF' and
% 'calRHVecBF' for the identification of body force.
fprintf(1,'Solving the linear system: \n');

startTime = cputime;

%fwdMapBFDir = [reslDir filesep 'fwdMapBF'];

%coefBFx = zeros(dimFS,numDTimePts);
%coefBFy = zeros(dimFS,numDTimePts);
ans = input('Select time steps (0 for all):');
if isempty(ans) || ans == 0
   selTimeSteps = 1:numDTimePts;
else
   selTimeSteps = ans;
end

for ii = 1:length(selTimeSteps)
   jj = selTimeSteps(ii);

   localStartTime = cputime;

   fprintf(1,'   Time step %d ...',jj);
   imgIndex = imgIndexOfDTimePts(jj);

   fwdMapBFFile = [fwdMapBFDir filesep 'A' sprintf(imgIndexForm,imgIndex) '.mat'];
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

   regParamFile = [reslDir filesep 'lastSavedRegParam.mat'];
   if exist(regParamFile,'file') == 2
      fprintf(1,'\n');
      s = load(regParamFile);
      regParam = s.regParam;

      answer = input(sprintf(['   Last saved optimal range of regularization parameter: ' ...
         '(%3.2f,%3.2f)\n' '   Do you want to reidentify it using L-curve? (y/n):'], ...
         regParam.minBFSigma,regParam.maxBFSigma),'s');
   else
      answer = 'yes';
   end

   %Identify the optimal range of regularization parameter using L-curve.
   if strcmp(answer,'n')
      answer = input(sprintf('   Choose your sigma in this range (last saved: %3.2f):', ...
         regParam.selBFSigma),'s');
      if strcmp(answer,'')
         selBFSigma = regParam.selBFSigma;
      else
         selBFSigma = str2num(answer);
         regParam.selBFSigma = selBFSigma;
         save(regParamFile,'regParam');
      end
      isBFSigmaRangeIdentified = 'yes';
   else
      fprintf(1,'   Identifying optimal sigma range using L-curve ...\n');
      isBFSigmaRangeIdentified = 'no';
   end

   while strcmp(isBFSigmaRangeIdentified,'no')
      %Calculate the L-curve.
      bfSigmaCandidate = bfSigma*[1e-4 1e-3 1e-2 1e-1 1 1e1 1e2 1e3 1e4];
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
         regParam.selBFSigma = str2num(answer);
         regParam.minBFSigma = regParam.selBFSigma*1e-1;
         regParam.maxBFSigma = regParam.selBFSigma*1e1;
         isBFSigmaRangeIdentified = 'yes';
         save(regParamFile,'regParam');
      else
         answer = input('   Please select a new center of test range:','s');
         bfSigma = str2num(answer);
         isBFSigmaRangeIdentified = 'no';
      end
   end

   B = B0 + regParam.selBFSigma*eye(n);
   coef = B\b;
   %coef = (A{jj}.'*A{jj}+bfSigma*eye(n))\(A{jj}.'*rightU{jj});

   forceField.coefBF = zeros(dimFS,2);
   forceField.coefBF(indDomDOF,1) = coef(1:dimBF);
   forceField.coefBF(indDomDOF,2) = coef(dimBF+1:2*dimBF);

   %Save the coefficient.
   forceFieldFile = [forceFieldDir filesep 'forceField' ...
      sprintf(imgIndexForm,imgIndex) '.mat'];
   save(forceFieldFile,'forceField');

   fprintf(1,'   Done in %5.3f sec.\n',cputime-localStartTime);
end

%Save the coefficient.
%save([reslDir filesep 'coefBFId'],'fs','coefBFx','coefBFy');

fprintf(1,'Total time spent: %5.3f sec.\n',cputime-startTime);

