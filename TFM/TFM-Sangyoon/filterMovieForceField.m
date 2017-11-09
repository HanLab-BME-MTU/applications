function [] = filterMovieForceField(MD)
% function [] = filterMovieForceField(MD) filters forceField for short
% peaks by applying filterForceShortPeaks which uses median filter.
% input: MD:    movieData file
% you have to have your TFM package run.

%% Input
% Get whole frame number
nFrames = MD.nFrames_;
% Get TFM package
TFMPackage = MD.getPackage(MD.getPackageIndex('TFMPackage'));
% Load Cell Segmentation
iMask = MD.getProcessIndex('MaskRefinementProcess');
if isempty(iMask)
    iMask = MD.getProcessIndex('ThresholdProcess');
end
if ~isempty(iMask)
    maskProc = MD.getProcess(iMask);
else
    error('You have to run ThresholdProcess or MaskRefinementProcess before running this code!')
end
% Load the forcefield
iForceFieldProc = 4;
forceFieldProc=TFMPackage.processes_{iForceFieldProc};
forceField=forceFieldProc.loadChannelOutput;
% save forceField
display('Backing up existing force field...')
tic
[pathstr1,name1,ext] = fileparts(forceFieldProc.outFilePaths_{1});
backUpPathTractionMap =  [pathstr1 filesep name1,'BeforeFiltered' ext];
if ~exist(backUpPathTractionMap,'file') && exist(forceFieldProc.outFilePaths_{1},'file')
    copyfile(forceFieldProc.outFilePaths_{1},backUpPathTractionMap)
end
toc
forceField = filterForceShortPeaks(forceField);

save(forceFieldProc.outFilePaths_{1},'forceField');

display('Backing up existing traction map...')
tic
[pathstr2,name2,ext] = fileparts(forceFieldProc.outFilePaths_{2});
backUpPathTractionMap =  [pathstr2 filesep name2,'BeforeFiltered' ext];
if ~exist(backUpPathTractionMap,'file') && exist(forceFieldProc.outFilePaths_{2},'file')
    copyfile(forceFieldProc.outFilePaths_{2},backUpPathTractionMap)
end
toc

iDisplFieldProc =MD.getProcessIndex('DisplacementFieldCorrectionProcess',1,0);     
if isempty(iDisplFieldProc)
    iDisplFieldProc =MD.getProcessIndex('DisplacementFieldCalculationProcess',1,0);     
end
displFieldProc=MD.processes_{iDisplFieldProc};
displField=displFieldProc.loadChannelOutput;

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
iBeadChan=1;
progressText(0,'Creating new traction map')
% Use mask of first frame to filter displacementfield
iSDCProc =MD.getProcessIndex('StageDriftCorrectionProcess',1,1);     
if ~isempty(iSDCProc)
    SDCProc=MD.processes_{iSDCProc};
    if ~SDCProc.checkChannelOutput(iBeadChan)
        error(['The channel must have been corrected ! ' ...
            'Please apply stage drift correction to all needed channels before '...
            'running displacement field calclation tracking!'])
    end
    %Parse input, store in parameter structure
    refFrame = double(imread(SDCProc.outFilePaths_{2,iBeadChan}));
    firstMask = false(size(refFrame));
else
    iDisplFieldProc =MD.getProcessIndex('DisplacementFieldCalculationProcess',1,0);     
    displFieldProc=MD.processes_{iDisplFieldProc};
    refFrame = imread(displFieldProc.funParams_.referenceFramePath);
    firstMask = false(size(refFrame));
end
for ii=1:nFrames
    % starts with original size of beads
    cur_tMap = zeros(size(firstMask));
    cur_tMapX = zeros(size(firstMask));
    cur_tMapY = zeros(size(firstMask));
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
    save(forceFieldProc.outFilePaths_{2},'tMap','tMapX','tMapY','fCfdMap'); % need to be updated for faster loading. SH 20141106
catch
    save(forceFieldProc.outFilePaths_{2},'tMap','tMapX','tMapY'); % need to be updated for faster loading. SH 20141106
end
forceFieldProc.setTractionMapLimits([tmin tmax])

