function [imgStack, tMap, imgStack2] = getAnyStacks(MD)
% Load FA package
FAPack=MD.getPackage(MD.getPackageIndex('FocalAdhesionPackage'));
% Load classification process
classProc = FAPack.getProcess(8);
% p = classProc.funParams_;

% Load the selectedGroups.mat
iChan = find(classProc.checkChannelOutput);

tfmPack = MD.getPackage(MD.getPackageIndex('TFMPackage'));
nFrames = MD.nFrames_; tMap=[];

% tMap creation
iBeadChan = 1; % might need to be updated based on asking TFMPackage..
SDCProc_FA= FAPack.processes_{1};
if ~isempty(SDCProc_FA)
    s = load(SDCProc_FA.outFilePaths_{3,iBeadChan},'T');    
    T_FA = s.T;
else
    T_FA = zeros(nFrames,2);
end

SDCProc_TFM=tfmPack.processes_{1};
%iSDCProc =MD.getProcessIndex('StageDriftCorrectionProcess',1,1);     
if ~isempty(SDCProc_TFM)
    s = load(SDCProc_TFM.outFilePaths_{3,iBeadChan},'T');    
    T_TFM = s.T;
else
    T_TFM = zeros(nFrames,2);
end
T = -T_TFM + T_FA;

for ii=nFrames:-1:1
    cur_tMap=tfmPack.processes_{4}.loadChannelOutput(ii,'output','tMap');
    cur_T = T(ii,:);
    cur_tMap2 = imtranslate(cur_tMap, cur_T(2:-1:1));
    tMap(:,:,ii) = cur_tMap2;
end

% Other image maps
sdcProc = FAPack.processes_{1};
imgStack = sdcProc.loadOutStack(2);
imgStack2 = sdcProc.loadOutStack(3);