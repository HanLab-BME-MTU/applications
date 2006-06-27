oldRawDispFieldFile = [mechDir filesep 'oldRawDispField.mat'];

if exist(oldRawDispFieldFile,'file') ~= 2
   fprintf(1,'There is no old raw displacement field saved.');
   return;
end

%Create the directory for raw displacement field if it does not exist.
rawDispFieldDir = [reslDir filesep 'rawDispField'];
if ~isdir(rawDispFieldDir)
   success = mkdir(mechDir,'rawDispField');
   if ~success
      error('Trouble making directory for raw displacement field.');
   end
end

load(oldRawDispFieldFile);
rawDispField.p       = rawDataP{1};
rawDispField.v       = rawDispV{1}(:,4:-1:3) - rawDataP{1};
rawDispField.outlier = rawOutlier{1};
rawDispField.inlier  = rawInlier{1};

rawDispFieldFileName = ['rawDispField' sprintf(DTIndexForm,1) '.mat'];
rawDispFieldFile     = [rawDispFieldDir filesep rawDispFieldFileName];
save(rawDispFieldFile,'rawDispField');


%Convert the old interpolated displacement field.
oldDispFieldFile = [reslDir filesep 'oldDispField.mat'];
if exist(oldDispFieldFile,'file') ~= 2
   fprintf(1,'There is no old interpolated displacement field saved\n.');
   return;
end

load(oldDispFieldFile);
iDispField.p         = [iDataPx{1} iDataPy{1}];
iDispField.v         = [iDataU1{1} iDataU2{1}];
iDispField.sv        = [sDataU1{1} sDataU2{1}];
iDispField.corLen    = MFCorLen;
iDispField.sCorLen   = corLen;
iDispField.edgCorLen = edgCorLen;
iDispField.gridx     = gridx;
iDispField.gridy     = gridy;
iDispField.gridX     = gridX;
iDispField.gridY     = gridY;
iDispField.gridPx    = gridPx;
iDispField.gridPy    = gridPy;

iDispFieldFileName = ['iDispField' sprintf(DTIndexForm,1) '.mat'];
iDispFieldFile     = [reslDir filesep 'iDispField' filesep iDispFieldFileName];
save(iDispFieldFile,'iDispField');

