function [] = adjustRelativeForceField(MD)
% function [] = adjustRelativeForceField(MD) adjusts forceField calculated
% based on the first bead frame to be all inward - based on cell mask.
% input: MD:    movieData file
% you have to have your TFM package run (based on the first frame of the
% bead channel), and segmentation package.

%% Input
% Get whole frame number
nFrames = MD.nFrames_;
% Get TFM package
TFMPackage = MD.getPackage(MD.getPackageIndex('TFMPackage'));
% Load Cell Segmentation
iSegPack=MD.getPackageIndex('SegmentationPackage');
if ~isempty(iSegPack)
    SegPackage = MD.getPackage(iSegPack);
    maskProc = SegPackage.getProcess(2);
else
    iMask = MD.getProcessIndex('MaskRefinementProcess');
    if isempty(iMask)
        iMask = MD.getProcessIndex('ThresholdProcess');
    end
    if ~isempty(iMask)
        maskProc = MD.getProcess(iMask);
    else
        error('You have to run ThresholdProcess or MaskRefinementProcess before running this code!')
    end
end
% Load the forcefield
iForceFieldProc = 4;
forceFieldProc=TFMPackage.processes_{iForceFieldProc};
forceField=forceFieldProc.loadChannelOutput;

% Channel for the cell mask is 2
iChan = 2;
iBeadChan = 1;

iSDCProc =MD.getProcessIndex('StageDriftCorrectionProcess',1,1);     
if ~isempty(iSDCProc)
    SDCProc=MD.processes_{iSDCProc};
    s = load(SDCProc.outFilePaths_{3,iBeadChan},'T');    
    T = s.T;
    maxX = ceil(max(abs(T(:, 2))));
    maxY = ceil(max(abs(T(:, 1))));
else
    T = zeros(nFrames,2);
    maxX = ceil(max(abs(T(:, 2))));
    maxY = ceil(max(abs(T(:, 1))));
end
% make stack of transforms
% for k=1:nFrames
k=1;
    cellMaskOrg = maskProc.loadChannelOutput(iChan,k);
    Tr = affine2d([1 0 0; 0 1 0; fliplr(T(k, :)) 1]);
    % Apply subpixel-wise registration to original masks
    I = padarray(cellMaskOrg, [maxY, maxX]);
    cellMask = imwarp(I, Tr);
    % Find the innermost point 
    segCentroid = regionprops(cellMask,'Centroid');
    iMPoint = segCentroid.Centroid; % this might be improved with putting more weight on the pixels at the boundary.
% end
% backup forceField
[pathstr1,name1,ext] = fileparts(forceFieldProc.outFilePaths_{1});
save([pathstr1 filesep name1,'Org' ext],'forceField');
%% Go through each node and get the vector
firstForceField = forceField(1).pos;
numNodes =length(firstForceField);
splineParam = 0.01;
progressText(0,'Vector adjusting')
filterWindow=3;
tRange = 1:nFrames;
for ii=1:numNodes
    % get the vector for entire time series
    curVecX = arrayfun(@(x) x.vec(ii,1),forceField);
    curVecY = arrayfun(@(x) x.vec(ii,2),forceField);
    % Smooth these out first

    curVecX_med = medfilt1(curVecX,filterWindow);
    curVecY_med = medfilt1(curVecY,filterWindow);
    curVecX_spline= csaps(tRange,curVecX_med,splineParam);
    curVecX_spline_discretized=ppval(curVecX_spline,tRange);
    curVecY_spline= csaps(tRange,curVecY_med,splineParam);
    curVecY_spline_discretized=ppval(curVecY_spline,tRange);
    
    % Find the vector to the node from the innermost point
%     Tr = maketform('affine', [1 0 0; 0 1 0; fliplr(T(ii, :)) 1]);
%     figure, imshow(cellMask), hold on, plot(iMPoint(1),iMPoint(2),'r*')
    curPosX = forceField(1).pos(ii,1);
    curPosY = forceField(1).pos(ii,2);
    vecToCurPos = [curPosX-iMPoint(1) curPosY-iMPoint(2)];
    vecToCurPos = vecToCurPos/norm(vecToCurPos); % normalize it
    % Dot product
%     curVec = [curVecX_med' curVecY_med'];
    curVec = [curVecX_spline_discretized' curVecY_spline_discretized'];
    dotProd = curVec * vecToCurPos';
%     figure, plot(1:nFrames,dotProd)
    % get indices containing positive values
    idxPositiveProd = find(dotProd>0);
%     figure, plot(tRange, curVecX), hold on, plot(tRange, curVecX_med)
%     figure, plot(tRange, curVecX), hold on, plot(tRange, curVecX_spline_discretized)
    % Among those, get the frame number that contains largest norm
    if ~isempty(idxPositiveProd)
%         curMag = sqrt(curVec(:,1).^2+curVec(:,2).^2);
%         [~,iNegFrameRel] = max(curMag(idxPositiveProd));
%         iNegFrame = idxPositiveProd(iNegFrameRel);
%         shiftVector = curVec(iNegFrame,:);
    %     % compare with frame contiaing most negative prod value
        [~,iMaxProd] = max(dotProd(idxPositiveProd));
        iMaxProd = idxPositiveProd(iMaxProd);
        % Shift the entire time series to make the vector zero
        shiftVector = curVec(iMaxProd,:);
        shiftedVectors = curVec - [shiftVector(1)*ones(nFrames,1) shiftVector(2)*ones(nFrames,1)];
        for k=1:nFrames
            forceField(k).vec(ii,:) = [shiftedVectors(k,1), shiftedVectors(k,2)];
        end
    end
    progressText(ii/(numNodes));
end
% save forceField
save(forceFieldProc.outFilePaths_{1},'forceField');

%% make traction map
iDisplFieldProc =MD.getProcessIndex('DisplacementFieldCorrectionProcess',1,0);     
if isempty(iDisplFieldProc)
    iDisplFieldProc =MD.getProcessIndex('DisplacementFieldCalculationProcess',1,0);     
end
displFieldProc=MD.processes_{iDisplFieldProc};
displField=displFieldProc.loadChannelOutput;
display('Backing up existing traction map...')
tic
[pathstr2,name2,ext] = fileparts(forceFieldProc.outFilePaths_{2});
backUpPathTractionMap =  [pathstr2 filesep name2,'Org' ext];
if ~exist(backUpPathTractionMap,'file') && exist(forceFieldProc.outFilePaths_{2},'file')
    copyfile(forceFieldProc.outFilePaths_{2},backUpPathTractionMap)
end
toc
display('Creating new traction map...')
tic
[tMapIn, tmax, tmin, cropInfo,tMapXin,tMapYin] = generateHeatmapShifted(forceField,displField,0);
display(['Estimated traction maximum = ' num2str(tmax) ' Pa.'])
forceFieldShifted(nFrames)=struct('pos','','vec','');
tMapAll=load(forceFieldProc.outFilePaths_{2});
tMap = tMapAll.tMap;
tMapX = tMapAll.tMapX;
tMapY = tMapAll.tMapY;
try
    fCfdMap = tMapAll.fCfdMap;
    buildForceConfidence=false;
catch
    try
        M=load(forceFieldProc.outFilePaths_{3},'M');
        M=M.M;
        forceNodeMaxima = max(M);
        forceConfidence.pos = forceMesh.p;
        forceConfidence.vec = reshape(forceNodeMaxima,[],2);
        %Make it relative
        maxCfd = max(forceNodeMaxima);
        forceConfidence.vec = forceConfidence.vec/maxCfd;
        [fCfdMapIn] = generateHeatmapShifted(forceConfidence,displField,0);
        fCfdMapIn{1} = fCfdMapIn{1}/max(fCfdMapIn{1}(:));
        buildForceConfidence=true;
    catch
        buildForceConfidence=false;
    end
end
toc
% distBeadMap = cell(1,nFrames);
%% assigning to the file
reg_grid=createRegGridFromDisplField(displField,1,0);
progressText(0,'Creating new traction map')
for ii=1:nFrames
    % starts with original size of beads
    cur_tMap = zeros(size(cellMask));
    cur_tMapX = zeros(size(cellMask));
    cur_tMapY = zeros(size(cellMask));
    cur_tMap(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = tMapIn{ii};
    cur_tMapX(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = tMapXin{ii};
    cur_tMapY(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = tMapYin{ii};
    tMap{ii} = cur_tMap;
    tMapX{ii} = cur_tMapX;
    tMapY{ii} = cur_tMapY;
    % Shifted forceField vector field
    curDispVec = displField(ii).vec;
    curDispPos = displField(ii).pos;
    curDispVec = curDispVec(~isnan(curDispVec(:,1)),:); % This will remove the warning 
    curDispPos = curDispPos(~isnan(curDispVec(:,1)),:); % This will remove the warning 
    [grid_mat,iu_mat, ~,~] = interp_vec2grid(curDispPos, curDispVec,[],reg_grid);
    displ_vec = [reshape(iu_mat(:,:,1),[],1) reshape(iu_mat(:,:,2),[],1)]; 
    
    [forceFieldShiftedpos,forceFieldShiftedvec, ~, ~] = interp_vec2grid(forceField(ii).pos+displ_vec, forceField(ii).vec,[],grid_mat); %1:cluster size
    pos = [reshape(forceFieldShiftedpos(:,:,1),[],1) reshape(forceFieldShiftedpos(:,:,2),[],1)]; %dense
    force_vec = [reshape(forceFieldShiftedvec(:,:,1),[],1) reshape(forceFieldShiftedvec(:,:,2),[],1)]; 
    if ii==1 && buildForceConfidence
        cur_fCfdMap = zeros(size(cellMask));
        cur_fCfdMap(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = fCfdMapIn{ii};
        fCfdMap = cur_fCfdMap;
    end     

    forceFieldShifted(ii).pos = pos;
    forceFieldShifted(ii).vec = force_vec;
    progressText(ii/(nFrames));
end
save(forceFieldProc.outFilePaths_{1},'forceField','forceFieldShifted');
try
    save(forceFieldProc.outFilePaths_{2},'tMap','tMapX','tMapY','fCfdMap','-v7.3'); % need to be updated for faster loading. SH 20141106
catch
    save(forceFieldProc.outFilePaths_{2},'tMap','tMapX','tMapY','-v7.3'); % need to be updated for faster loading. SH 20141106
end
forceFieldProc.setTractionMapLimits([tmin tmax])

