function [imgStack, tMap, imgStack2] = getAnyStacks(MD)
% Load FA package
FAPack=MD.getPackage(MD.getPackageIndex('FocalAdhesionPackage'));
% Load classification process
classProc = FAPack.getProcess(8);
% p = classProc.funParams_;

% Load the selectedGroups.mat
iChan = find(classProc.checkChannelOutput);

iTFMPackage=MD.getPackageIndex('TFMPackage');
if ~isempty(iTFMPackage)
    tfmPack = MD.getPackage(iTFMPackage);
end
nFrames = MD.nFrames_; tMap=[];

% tMap creation
iBeadChan = 1; % might need to be updated based on asking TFMPackage..
SDCProc_FA= FAPack.processes_{1};
if ~isempty(SDCProc_FA)
    try
        s = load(SDCProc_FA.outFilePaths_{3,iBeadChan},'T');    
        T_FA = s.T;
    catch
        T_FA = zeros(nFrames,2);
    end
else
    T_FA = zeros(nFrames,2);
end

if ~isempty(iTFMPackage)
    SDCProc_TFM=tfmPack.processes_{1};
    %iSDCProc =MD.getProcessIndex('StageDriftCorrectionProcess',1,1);     
    if ~isempty(SDCProc_TFM)
        s = load(SDCProc_TFM.outFilePaths_{3,iBeadChan},'T');    
        T_TFM = s.T;
    else
        T_TFM = zeros(nFrames,2);
    end
else
    T_TFM = zeros(nFrames,2);
end

T = -T_TFM + T_FA;

if ~isempty(iTFMPackage)
    if nargout>1 && ~isempty(tfmPack.processes_{4})
        tMapOrg=tfmPack.processes_{4}.loadChannelOutput('output','tMap');
        for ii=nFrames:-1:1
            cur_tMap=tMapOrg(:,:,ii);
            cur_T = T(ii,:);
            cur_tMap2 = imtranslate(cur_tMap, cur_T(2:-1:1));
            tMap(:,:,ii) = cur_tMap2;
        end
    end
end

% Other image maps
if ~isempty(SDCProc_FA)
    if ismember(2, find(SDCProc_FA.checkChannelOutput))
        imgStack = SDCProc_FA.loadOutStack(iChan);
    else
        imgStack = [];
    end
else
    imgStack = MD.channels_(iChan).loadImage(1:nFrames);
end

if nargout>2
    if ~isempty(SDCProc_FA)
        if ismember(iChan+1, find(SDCProc_FA.checkChannelOutput)) && ...
                ~isempty(FAPack.getProcess(9))
            imgStack2 = SDCProc_FA.loadOutStack(iChan+1);
        else
            imgStack2 = [];
        end
    else
        if numel(MD.channels_)>2
            imgStack2 = MD.channels_(iChan+1).loadImage(1:nFrames);
        else
            imgStack2 = [];
        end
    end
end