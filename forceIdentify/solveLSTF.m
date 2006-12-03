%This script file solves the linear system assembled by 'calFwdOpTF' and
% 'calRHVecTF' for the identification of the boundary traction force.
fprintf(1,'Solving the linear system: \n');

regParamFile = [reslDir filesep 'lastSavedRegParam.mat'];
if exist(regParamFile,'file') ~= 2
   fprintf(1,['Regularization parameter has not been identified yet.\n' ...
      'Please run ''calOptRegParTF'' first.\n']);
   return;
end

s = load(regParamFile);
regParam = s.regParam;

if ~isfield(regParam,'selTFSigma')
   fprintf(1,'Regularization parameter has not been identified yet.\n',regParam.selTFSigma);
   return;
end

startTime = cputime;

%Regularization parameter for 'solveLS'.
%sigma = 1;

answer = input('Select time steps (0 for all):');
if isempty(answer) | answer == 0
   selTimeSteps = 1:numDTimePts;
else
   selTimeSteps = answer;
end

for ii = 1:length(selTimeSteps)
   clear A;
   jj = selTimeSteps(ii);
   imgIndex = imgIndexOfDTimePts(jj);

   localStartTime = cputime;

   fprintf(1,'   Time step %d ...',jj);

   clear forceField;

   forceFieldFile = [forceFieldDir filesep 'forceField' ...
      sprintf(imgIndexForm,imgIndex) '.mat'];

   selTFSigma = [];
   if exist(forceFieldFile,'file') == 2
      s = load(forceFieldFile);
      forceField = s.forceField;

      if isfield(forceField,'selTFSigma')
         bfSigma    = forceField.bfSigma;
         selTFSigma = forceField.selTFSigma;
         minTFSigma = forceField.minTFSigma;
         maxTFSigma = forceField.maxTFSigma;

         regParam.bfSigma(jj)    = bfSigma;
         regParam.selTFSigma(jj) = selTFSigma;
         regParam.minTFSigma(jj) = minTFSigma;
         regParam.maxTFSigma(jj) = maxTFSigma;

         save(regParamFile,'regParam');
      end
   end

   if isempty(selTFSigma) && jj > 1
      %Try load the previous time point and use its regularization parameter.
      preImgIndex = imgIndexOfDTimePts(jj-1);
      preForceFieldFile = [forceFieldDir filesep 'forceField' ...
         sprintf(imgIndexForm,preImgIndex) '.mat'];
      if exist(preForceFieldFile,'file') == 2
         s = load(preForceFieldFile);
         preForceField = s.forceField;

         if isfield(preForceField,'selTFSigma')
            bfSigma    = preForceField.bfSigma;
            selTFSigma = preForceField.selTFSigma;
            minTFSigma = preForceField.minTFSigma;
            maxTFSigma = preForceField.maxTFSigma;

            forceField.bfSigma    = bfSigma;
            forceField.selTFSigma = selTFSigma;
            forceField.minTFSigma = minTFSigma;
            forceField.maxTFSigma = maxTFSigma;

            regParam.bfSigma(jj)    = bfSigma;
            regParam.selTFSigma(jj) = selTFSigma;
            regParam.minTFSigma(jj) = minTFSigma;
            regParam.maxTFSigma(jj) = maxTFSigma;

            save(forceFieldFile,'forceField');
            save(regParamFile,'regParam');
         end
      end

   end

   if isempty(selTFSigma)
      fprintf(1,['   Regularization parameter has not been identified yet.\n' ...
         '   Please run ''calOptRegParTF'' first.\n']);
      return;
   end


   fprintf(1,'   Regularization parameter: %5.3f.\n',selTFSigma);

   if strcmp(isFieldBndFixed,'yes') && strcmp(isDataPosFixed,'yes')
      fwdMapTFFile = [fwdMapTFDir filesep 'A' sprintf(imgIndexForm,0) '.mat'];
   else
      fwdMapTFFile = [fwdMapTFDir filesep 'A' sprintf(imgIndexForm,imgIndex) '.mat'];
   end
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
   numEdges = femModel.numEdges;

   procStr = '';
   for k = 1:numEdges
      %Reshape 'A';
      A{k} = reshape(A{k},2*numDP,2*fsBnd(k).dim);
      [m,n] = size(A{k});

      forceField.coefTF{k} = (A{k}.'*A{k}+selTFSigma*eye(n))\(A{k}.'*rightU{k});
      for kk = 1:length(procStr)
         fprintf(1,'\b');
      end
      procStr = sprintf('   Edge %d (out of %d) finished in %5.3f sec.', ...
         k,numEdges,cputime-localStartTime);
      fprintf(1,procStr);
   end
   fprintf(1,'\n');

   save(forceFieldFile,'forceField');
end

fprintf(1,'Total time spent: %5.3f sec.\n',cputime-startTime);

