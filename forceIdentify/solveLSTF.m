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

fprintf(1,'   Regularization parameter: %5.3f.\n',regParam.selTFSigma);
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

   localStartTime = cputime;

   fprintf(1,'   Time step %d ...',jj);
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

   forceFieldFile = [forceFieldDir filesep 'forceField' ...
      sprintf(imgIndexForm,imgIndex) '.mat'];
   s = load(forceFieldFile);
   forceField = s.forceField;

   procStr = '';
   for k = 1:numEdges
      %Reshape 'A';
      A{k} = reshape(A{k},2*numDP,2*fsBnd(k).dim);
      [m,n] = size(A{k});

      forceField.coefTF{k} = (A{k}.'*A{k}+regParam.selTFSigma*eye(n))\(A{k}.'*rightU{k});
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

