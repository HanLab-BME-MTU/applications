function [fsmParam,status,speckleArray]=fsmBuildMain(fsmParam)
% fsmBuildMain is the main function of the fsmBuild module
%
% SYNOPSIS   [fsmParam,status,speckleArray]=fsmBuildMain(fsmParam)
%
% INPUT      fsmParam     :   general parameter structure
%
% OUTPUT     speckleArray :   structure containing all speckle information from a movie
%            fsmParam     :   modified (when needed) parameter structure
%            status       :   it decides whether the program should continue after this module 
%                             or should stop because of errors;
%                             status is set to 0 (error) in the beginning of a module and will
%                             be set to 1 at the end if the module completed successfully.
%                             status = 1 - if the module completed successfully
%                             status = 0 - if the module did not complete successfully
%
% DEPENDENCES   fsmBuildMain uses {}
%               fsmBuildMain is used by { fsmMain }
%
% Aaron Ponti, October 7th, 2002

% Set initial module status
status=0;

% Check input parameter
if nargin~=1
    error('Input parameter fsmParam expected');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% READ NEEDED PARAMETERS FROM fsmParam
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xmin=fsmParam.main.normMin;
xmax=fsmParam.main.normMax;
threshold=fsmParam.track.threshold;
noiseParam=fsmParam.main.noiseParam;
userPath=fsmParam.main.path;
outFileList=fsmParam.specific.fileList;
strg=fsmParam.specific.formString;
factors=fsmParam.specific.intCorrFactors;
imageSize=fsmParam.specific.imgSize;
W=[1 1 0 0; 0 0 imageSize(1) imageSize(2)]; % Maintained for possible future re-use
firstIndex=fsmParam.specific.firstIndex;
sigma=fsmParam.prep.sigma;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CHANGE TO WORK PATH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(userPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% LOAD MAGIC POSITION MATRIX FROM DISK
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist([userPath,filesep,'mpm.mat'])==2
    load([userPath,filesep,'mpm.mat']);
else
    error('The file ''mpm.mat'' cannot be found. Please re-run the TRACKER module');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BUILD speckleArray
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save information into an array of struct speckle
% To cope with subpixel.
if fsmParam.prep.subpixel
    speckleArray=fsmBuildSaveSpeckleArraySPix(MPM,outFileList,xmin,xmax,firstIndex,strg,noiseParam,threshold,factors,W,userPath,sigma);
else
    speckleArray=fsmBuildSaveSpeckleArray(MPM,outFileList,xmin,xmax,firstIndex,strg,noiseParam,threshold,factors,W,userPath,sigma);
end
%Force it to 'speckleArray.spPos' to be double.
%speckleArray=fsmBuildSaveSpeckleArraySPix(MPM,outFileList,xmin,xmax,firstIndex,strg,noiseParam,threshold,factors,W,userPath,sigma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   SAVE speckleArray STRUCTURE
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

infoH=fsmGuiInfo; ch=get(infoH,'Children');
if length(userPath)>40
    strUserPath=[userPath(1:37),'...'];
else
    strUserPath=userPath;
end
textString=['Saving speckleArray to ',strUserPath,'\speckleArray.mat'];
set(ch(3),'String','BUILDER MODULE');
set(ch(2),'String',textString);
set(ch(1),'String','Prease wait...');
save speckleArray.mat speckleArray;
set(ch(1),'String','Done!');
pause(1);
close(infoH);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SETTING MODULE STATUS TO 1 AND RETURNING
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the status to 1 to mean that the module successfully finished
status=1;
