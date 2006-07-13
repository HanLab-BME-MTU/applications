%Identify the optimal regularization parameter range and choose the regularization parameter in this
%range. We use time step 1.
imgIndex = imgIndexOfDTimePts(1);

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
   answer = input(sprintf('   Choose your regularization parameter in this range (last saved: %3.2f):', ...
   regParam.selBFSigma),'s');
   if strcmp(answer,'')
      selBFSigma = regParam.selBFSigma;
   else
      selBFSigma = str2num(answer);
      regParam.selBFSigma = selBFSigma;
      save(regParamFile,'regParam');
   end
   isBFSigmaRangeIdentified = 'yes';
   return;
end

fprintf(1,'   Identifying optimal regularization parameter range using L-curve ...\n');
isBFSigmaRangeIdentified = 'no';

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

fprintf(1,'Optimal regularization parameter range: (%3.2f,%3.2f).\n', ...
   regParam.minBFSigma,regParam.maxBFSigma);
answer = input(sprintf('Choose your regularization parameter in this range (last saved: %3.2f):', ...
   regParam.selBFSigma),'s');
if strcmp(answer,'')
   selBFSigma = regParam.selBFSigma;
else
   selBFSigma = str2num(answer);
   regParam.selBFSigma = selBFSigma;
   save(regParamFile,'regParam');
end
isBFSigmaRangeIdentified = 'yes';

