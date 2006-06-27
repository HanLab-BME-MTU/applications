function setupForceProj(varargin)
%This script file set up the project for force reconstruction. It should be
% called first.
%
% SYNOPSIS:
%    setupForceProj;
%    setupForceProj('reslDir',reslDir,'simuDir',simuDir);
% INPUT:
%    Optional PAR/VALUE pairs:
%       PAR         VALUE
%    ---------------------------------------------------------------------
%    'reslDir' : The directory for saveing results. It should start with 'resl'.
%                By default, it is 'resl'.
%    'simuDir' : The directory for simulating retrograde flow. It should start 
%                with 'simu'. By default, it is [].

[projDir,imgDir,subProjDir,imgDirList,firstImgList] = ...
   projSetupGUI([],[]);

if isempty(projDir)
   fprintf('Project setup is canceled.\n');
   return;
end

%Default parameters.
reslDir = 'resl';
simuDir = [];

if mod(nargin,2) == 1
   error('Unbalanced PAR/VALUE pairs.');
end

for k = 1:2:nargin
   switch varargin{k}
      case 'reslDir'
         reslDir = varargin{k+1};
      case 'simuDir'
         simuDir = varargin{k+1};
   end
end

numImgChannels = length(firstImgList);
for k = 1:numImgChannels
    firstImgFileName = [imgDirList{k} filesep firstImgList{k}];
    imgFileList{k} = getFileStackNames(firstImgFileName);
end

%Get the first and last image index.
[path,body,no,ext] = getFilenameBody(imgFileList{1}{1});
firstImgIndex = str2num(no);

[path,body,no,ext] = getFilenameBody(imgFileList{1}{end});
lastImgIndex = str2num(no);

%Get the number (or index) format of the last image.
imgIndexForm  = sprintf('%%.%dd',length(num2str(no)));

tackDir = subProjDir{strmatch('tack',subProjDir)};
edgeDir = subProjDir{strmatch('edge',subProjDir)};
corrDir = subProjDir{strmatch('corr',subProjDir)};
mechDir = subProjDir{strmatch('mech',subProjDir)};

tackDir = [projDir filesep tackDir];
edgeDir = [projDir filesep edgeDir];
corrDir = [projDir filesep corrDir];
mechDir = [projDir filesep mechDir];

if ~isdir([mechDir filesep reslDir])
   mkdir(mechDir,reslDir);
end
reslDir = [mechDir filesep reslDir];

if ~isempty(simuDir)
   if ~isdir([mechDir filesep simuDir])
      mkdir(mechDir,simuDir);
   end
   simuDir = [mechDir filesep simuDir];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create all subdirectories under 'mechDir' or 'reslDir' for storing results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create the directory for raw displacement field if it does not exist.
rawDispFieldDir = [reslDir filesep 'rawDispField'];
if ~isdir(rawDispFieldDir)
   success = mkdir(reslDir,'rawDispField');
   if ~success
      error('Trouble making directory for raw displacement field.');
   end
end

%Create the directory for storing field boundary if it does not exist.
%fieldBndDir = [mechDir filesep 'fieldBnd'];
%if ~isdir(fieldBndDir)
%   success = mkdir(mechDir,'fieldBnd');
%   if ~success
%      error('Trouble making directory for field boundary.');
%   end
%end

%Create the directory for interpolated displacement field if it does not exist.
iDispFieldDir = [reslDir filesep 'iDispField'];
if ~isdir(iDispFieldDir)
   success = mkdir(reslDir,'iDispField');
   if ~success
      error('Trouble making directory for interpolated displacement field.');
   end
end

%Create the directory for the fem model if it does not exist.
femModelDir = [mechDir filesep 'femModel'];
if ~isdir(femModelDir)
   success = mkdir(mechDir,'femModel');
   if ~success
      error('Trouble making directory for the fem model.');
   end
end

femSolBasisBFDir = [mechDir filesep 'femSolBasisBF'];
if ~isdir(femSolBasisBFDir)
   success = mkdir(mechDir,'femSolBasisBF');
   if ~success
      error('Trouble making directory for ''femSolBasisBF''.');
   end
end

fwdMapBFDir = [reslDir filesep 'fwdMapBF'];
if ~isdir(fwdMapBFDir)
   success = mkdir(reslDir, 'fwdMapBF');
   if ~success
      error('Trouble making directory for ''fwdMapBF''.');
   end
end

femSolBasisTFDir = [mechDir filesep 'femSolBasisTF'];
if ~isdir(femSolBasisTFDir)
   success = mkdir(mechDir,'femSolBasisTF');
   if ~success
      error('Trouble making directory for ''femSolBasisTF''.');
   end
end

fwdMapTFDir = [reslDir filesep 'fwdMapTF'];
if ~isdir(fwdMapTFDir)
   success = mkdir(reslDir, 'fwdMapTF');
   if ~success
      error('Trouble making directory for ''fwdMapTF''.');
   end
end

forceFieldDir = [reslDir filesep 'forceField'];
if ~isdir(forceFieldDir)
   success = mkdir(reslDir,'forceField');
   if ~success
      error('Trouble making directory for interpolated displacement field.');
   end
end


%Check whether the two parameter files exist.
if exist([mechDir filesep 'modelPar.m'],'file') ~= 2
   defModelParFile = which('defModelPar.m');
   copyfile(defModelParFile,[mechDir filesep 'modelPar.m']);
end

if exist([reslDir filesep 'setPar.m'],'file') ~= 2
   setDefParFile = which('setDefPar.m');
   copyfile(setDefParFile,[reslDir filesep 'setPar.m']);
end

%Set the parameters.
resetPar;

