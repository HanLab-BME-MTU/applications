%This script file solves the linear system assembled by 'calFwdOpBF' and
% 'calRHVecBF' for the identification of body force.
fprintf(1,'Solving the linear system: \n');

regParamFile = [reslDir filesep 'lastSavedRegParam.mat'];
if exist(regParamFile,'file') ~= 2
   fprintf(1,['Regularization parameter has not been identified yet.\n' ...
      'Please run ''calOptRegParBF'' first.\n']);
   return;
end

s = load(regParamFile);
regParam = s.regParam;

startTime = cputime;

%fwdMapBFDir = [reslDir filesep 'fwdMapBF'];

%coefBFx = zeros(dimFS,numDTimePts);
%coefBFy = zeros(dimFS,numDTimePts);
answer = input('Select time steps (0 for all):');
if isempty(answer) | answer == 0
   selTimeSteps = 1:numDTimePts;
else
   selTimeSteps = answer;
end

for ii = 1:length(selTimeSteps)
   jj = selTimeSteps(ii);
   imgIndex = imgIndexOfDTimePts(jj);

   localStartTime = cputime;

   fprintf(1,'Time step %d ...',jj);

   clear forceField;

   forceFieldFile = [forceFieldDir filesep 'forceField' ...
      sprintf(imgIndexForm,imgIndex) '.mat'];

   selBFSigma = [];
   if exist(forceFieldFile,'file') == 2
      s = load(forceFieldFile);
      forceField = s.forceField;

      if isfield(forceField,'selBFSigma')
         bfSigma    = forceField.bfSigma;
         selBFSigma = forceField.selBFSigma;
         minBFSigma = forceField.minBFSigma;
         maxBFSigma = forceField.maxBFSigma;

         regParam.bfSigma(jj)    = bfSigma;
         regParam.selBFSigma(jj) = selBFSigma;
         regParam.minBFSigma(jj) = minBFSigma;
         regParam.maxBFSigma(jj) = maxBFSigma;

         save(regParamFile,'regParam');
      end
   end

   if isempty(selBFSigma) && jj > 1
      %Try load the previous time point and use its regularization parameter.
      preImgIndex = imgIndexOfDTimePts(jj-1);
      preForceFieldFile = [forceFieldDir filesep 'forceField' ...
         sprintf(imgIndexForm,preImgIndex) '.mat'];
      if exist(preForceFieldFile,'file') == 2
         s = load(preForceFieldFile);
         preForceField = s.forceField;

         if isfield(preForceField,'selBFSigma')
            bfSigma    = preForceField.bfSigma;
            selBFSigma = preForceField.selBFSigma;
            minBFSigma = preForceField.minBFSigma;
            maxBFSigma = preForceField.maxBFSigma;

            forceField.bfSigma    = bfSigma;
            forceField.selBFSigma = selBFSigma;
            forceField.minBFSigma = minBFSigma;
            forceField.maxBFSigma = maxBFSigma;

            regParam.bfSigma(jj)    = bfSigma;
            regParam.selBFSigma(jj) = selBFSigma;
            regParam.minBFSigma(jj) = minBFSigma;
            regParam.maxBFSigma(jj) = maxBFSigma;

            save(forceFieldFile,'forceField');
            save(regParamFile,'regParam');
         end
      end

   end

   if isempty(selBFSigma)
      fprintf(1,['   Regularization parameter has not been identified yet.\n' ...
         '   Please run ''calOptRegParBF'' first.\n']);
      return;
   end

   fprintf(1,'   Regularization parameter: %5.3f.\n',selBFSigma);

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

   B = B0 + selBFSigma*eye(n);
   coef = B\b;
   %coef = (A{jj}.'*A{jj}+bfSigma*eye(n))\(A{jj}.'*rightU{jj});

   forceField.coefBF = zeros(dimFS,2);
   forceField.coefBF(indDomDOF,1) = coef(1:dimBF);
   forceField.coefBF(indDomDOF,2) = coef(dimBF+1:2*dimBF);

   save(forceFieldFile,'forceField');

   fprintf(1,'   Done in %5.3f sec.\n',cputime-localStartTime);
end

%Save the coefficient.
%save([reslDir filesep 'coefBFId'],'fs','coefBFx','coefBFy');

fprintf(1,'Total time spent: %5.3f sec.\n',cputime-startTime);

