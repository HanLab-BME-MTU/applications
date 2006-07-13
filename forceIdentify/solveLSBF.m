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

fprintf(1,'   Regularization parameter: %5.3f.\n',regParam.selBFSigma);
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

