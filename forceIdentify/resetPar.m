%Rerun 'setPar' to reset the parameters. It can only run after you
% 'setupForceProj' so that the 'resl' directory is known and it knows where 
% to look for 'setPar'

run([reslDir filesep 'setPar.m']);

switch trackMethod
   case 'speck'
      dataFile = [tackDir filesep 'mpm.mat'];

      startImgIndex = firstImgIndex;
      endImgIndex   = startImgIndex+20;
      timeStepSize  = 0;
      numAvgFrames  = 20;
      curDTimePt    = 1;

      if timeStepSize == 0
         imgIndexOfDTimePts = startImgIndex;
      else
         imgIndexOfDTimePts = ...
            startImgIndex:timeStepSize:endImgIndex-numAvgFrames+1;
      end
   case 'dArray'
      curDTimePt = 1;
      dArrayDir = [tackDir filesep 'dArray'];
      [dispArrayFiles imgIndexOfDTimePts] = getDispArrayFiles(dArrayDir);
      startImgIndex = imgIndexOfDTimePts(1);
      endImgIndex   = imgIndexOfDTimePts(end);
      numAvgFrames  = ones(size(imgIndexOfDTimePts));

      if isempty(imgIndexOfDTimePts)
         fprintf(1,'No ''disp_array'' data is available.\n');
         dataFile = [];
      else
         for k = 1:length(dispArrayFiles)
            dataFile{k} = [dArrayDir filesep dispArrayFiles{k}];
         end
      end
   case 'corr'
      curDTimePt = 1;

      %Get all the flow track files.
      [flowTrackFile,imgIndexOfDTimePts,corrEndImgIndex] = ...
         getFlowTrackFiles(corrDir);

      if isempty(imgIndexOfDTimePts)
         fprintf(1,'No ''flowTrack'' data is available.\n');
         dataFile = [];
      else
         %If there are multiple 'flowTrack' data for one time point, only keep
         % the first one.
         iDiff = diff(imgIndexOfDTimePts);
         zeroI = find(iDiff==0);
         if ~isempty(zeroI)
            msg = sprintf(['Multiple ''flowTrack'' data files are ' ...
            'found for the following image index: %d ' ...
            'The first one is used.'], zeroI);
            warndlg(msg,'flowTrack data file','modal');

            imgIndexOfDTimePts(zeroI+1) = [];
            corrEndImgIndex(zeroI+1)    = [];
            flowTrackFile(zeroI+1)      = [];
         end

         for k = 1:length(flowTrackFile)
            dataFile{k} = [corrDir filesep flowTrackFile{k}];
         end
         startImgIndex = imgIndexOfDTimePts(1); 
         endImgIndex   = corrEndImgIndex(end); 
         numAvgFrames  = corrEndImgIndex-imgIndexOfDTimePts+1;
      end
end

if ~isempty(dataFile)
   numDTimePts = length(imgIndexOfDTimePts);
   DTIndexForm = sprintf('%%.%dd',max(2,length(num2str(numDTimePts))));

   %Get the first image of the cell.
   relFrameNo = imgIndexOfDTimePts(1)-firstImgIndex+1;
   firstDTImg = imread(imgFileList{1}{relFrameNo});
else
   numDTimePts = 0;
   DTIndexForm = '%.2d';
   relFrameNo = 1;
   firstDTImg = [];
   fprintf(1,'There is no displacement data file.\n');
end

%run([mechDir filesep 'modelPar.m']);

%Export all the variables.
allVar = who;
for k = 1:length(allVar)
   assignin('base',allVar{k},eval(allVar{k}));
end

%Save parameters
save([reslDir filesep 'lastSavedParam.mat'],'param');

%Save the project setting in 'resl' directoty.
load([projDir filesep 'lastProjSettings.mat']);

projSettings.dataFile           = dataFile;
projSettings.imgIndexOfDTimePts = imgIndexOfDTimePts;
projSettings.DTIndexForm        = DTIndexForm;
projSettings.imgIndexForm       = imgIndexForm;
projSettings.imgFileList        = imgFileList;

save([reslDir filesep 'lastProjSettings.mat'],'projSettings');

%Bundle all the parameters into 'fHandles'.
%fHandles.projDir = projDir;
%fHandles.imgDir  = imgDir;
%fHandles.tackDir = tackDir;
%fHandles.corrDir = corrDir;
%fHandles.mechDir = mechDir;
%fHandles.reslDir = reslDir;
%
%fHandles.trackMethod   = trackMethod;
%fHandles.dataFile      = dataFile;
%fHandles.indexForm     = indexForm;
%fHandles.firstImgIndex = firstImgIndex;
%fHandles.lastImgIndex  = lastImgIndex;
%fHandles.startImgIndex = startImgIndex;
%fHandles.endImgIndex   = endImgIndex;
%fHandles.numAvgFrames  = numAvgFrames;
%fHandles.imgFileList   = imgFileList;
%
%fHandles.imgIndexOfDTimePts = imgIndexOfDTimePts;
%fHandles.numDTimePts        = length(imgIndexOfDTimePts);
%fHandles.curDTimePt         = curDTimePt;
%
%fHandles.MFCorLen        = MFCorLen;
%fHandles.corLen          = corLen;
%fHandles.gridDx          = gridDx;
%fHandles.gridDy          = gridDy;
%fHandles.edgeCorLen      = edgeCorLen;
%fHandles.dispType        = dispType;
%fHandles.calInterp       = calInterp;
%fHandles.showInterp      = showInterp;
%fHandles.dataSite        = 'speckle;
%fHandles.edgeCorLen      = edgeCorLen;
%fHandles.fowdOpComputed  = fwdOpComputed;
%fHandles.forceToIdentify = forceToIdentify;
%fHandles.bfDisplaySite   = bfDispSite;
%fHandles.bfScale         = bfScale;
%fHandles.dispScale       = dispScale;
%fHandles.smDispThreshold = smDispThreshold;
%fHandles.mcfAngle        = mcfAngle;
%fHandles.adfAngle        = adfAngle;
%fHandles.dataToUse       = dataToUse;
%fHandles.edgeBrkDistTF   = edgeBrkDistTF;
%fHandles.bspOrderTF      = bspOrderTF;
%fHandles.sigma           = sigma;
%fHandles.testTimeStep    = testTimeStep;
%
