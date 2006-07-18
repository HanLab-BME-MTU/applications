%For boundary force:
%Identify the optimal regularization parameter range and choose the regularization parameter in this
%range. We use time step 1.
imgIndex = imgIndexOfDTimePts(1);

regParamFile = [reslDir filesep 'lastSavedRegParam.mat'];
if exist(regParamFile,'file') == 2
   fprintf(1,'\n');
   s = load(regParamFile);
   regParam = s.regParam;

   if isfield(regParam,'minTFSigma') && isfield(regParam,'maxTFSigma')
      answer = input(sprintf(['   Last saved optimal range of regularization parameter: ' ...
         '(%3.2f,%3.2f)\n' '   Do you want to reidentify it using L-curve? (y/n):'], ...
         regParam.minTFSigma,regParam.maxTFSigma),'s');
   else
      answer = 'yes';
   end
else
   answer = 'yes';
end

%Identify the optimal range of regularization parameter using L-curve.
if strcmp(answer,'n')
   answer = input(sprintf('   Choose your regularization parameter in this range (last saved: %3.2f):', ...
   regParam.selTFSigma),'s');
   if strcmp(answer,'')
      selTFSigma = regParam.selTFSigma;
   else
      selTFSigma = str2num(answer);
      regParam.selTFSigma = selTFSigma;
      save(regParamFile,'regParam');
   end
   isTFSigmaRangeIdentified = 'yes';
   return;
end

fprintf(1,'   Identifying optimal regularization parameter range using L-curve ...\n');
isTFSigmaRangeIdentified = 'no';

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
   tfSigmaCandidate = tfSigma*[1e-4 1e-3 1e-2 1e-1 1 1e1 1e2 1e3 1e4];

   coefNorm    = zeros(1,length(tfSigmaCandidate));
   residueNorm = zeros(1,length(tfSigmaCandidate));
   for jj = 1:numEdges
      %Reshape 'A';
      A{jj} = reshape(A{jj},2*numDP,2*fsBnd(jj).dim);
      [m,n] = size(A{jj});

      B0 = A{jj}.'*A{jj};
      b  = A{jj}.'*rightU{jj};

      %Calculate the L-curve.
      tfSigmaCandidate = tfSigma*[1e-4 1e-3 1e-2 1e-1 1 1e1 1e2 1e3 1e4];
      for kk = 1:length(tfSigmaCandidate)
         B_i = B0 + tfSigmaCandidate(kk)*eye(n);
         coef_i = B_i\b;
         coefNorm(kk)    = coefNorm(kk)+sqrt(sum(coef_i.^2));
         residueNorm(kk) = residueNorm(kk)+sqrt(sum((A{jj}*coef_i-rightU{jj}).^2));
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
      regParam.selTFSigma = str2num(answer);
      regParam.minTFSigma = regParam.selTFSigma*1e-1;
      regParam.maxTFSigma = regParam.selTFSigma*1e1;
      isTFSigmaRangeIdentified = 'yes';
      save(regParamFile,'regParam');
   else
      answer = input('   Please select a new center of test range:','s');
      tfSigma = str2num(answer);
      isTFSigmaRangeIdentified = 'no';
   end
end

fprintf(1,'Optimal regularization parameter range: (%3.2f,%3.2f).\n', ...
   regParam.minTFSigma,regParam.maxTFSigma);
answer = input(sprintf('Choose your regularization parameter in this range (last saved: %3.2f):', ...
   regParam.selTFSigma),'s');
if strcmp(answer,'')
   selTFSigma = regParam.selTFSigma;
else
   selTFSigma = str2num(answer);
   regParam.selTFSigma = selTFSigma;
   save(regParamFile,'regParam');
end
isTFSigmaRangeIdentified = 'yes';

