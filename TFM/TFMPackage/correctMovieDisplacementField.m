function correctMovieDisplacementField(movieData,varargin)
% correctMovieDisplacementField calculate the displacement field
%
% correctMovieDisplacementField 
%
% SYNOPSIS correctMovieDisplacementField(movieData,paramsIn)
%
% INPUT   
%   movieData - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%

% Sebastien Besson, Sep 2011
% Sangyoon Han, from Oct 2014
%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous stage drift processes                                                                     
iProc = movieData.getProcessIndex('DisplacementFieldCorrectionProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(DisplacementFieldCorrectionProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end
displFieldCorrProc = movieData.processes_{iProc};
%Parse input, store in parameter structure
p = parseProcessParams(displFieldCorrProc,paramsIn);

% Default values for 3D spatiotemporal filter parameters
% (added here for backward compatibility with older Process class definitions)
if ~isfield(p,'temporalWindow'),   p.temporalWindow   = 2;   end  % frames before/after
if ~isfield(p,'outlierThreshold3D'), p.outlierThreshold3D = 3; end  % MAD multiplier
if ~isfield(p,'spatialWeight3D'),  p.spatialWeight3D  = 0.6; end  % spatial vs temporal blend
if ~isfield(p,'temporalWeight3D'), p.temporalWeight3D = 0.4; end

%% Backup the original vectors to backup folder
if exist(p.OutputDirectory,'dir')
    disp('Backing up the original data')
    ii = 1;
    backupFolder = [p.OutputDirectory ' Backup ' num2str(ii)];
    while exist(backupFolder,'dir')
        backupFolder = [p.OutputDirectory ' Backup ' num2str(ii)];
        ii=ii+1;
    end
    mkdir(backupFolder);
    copyfile(p.OutputDirectory, backupFolder,'f')
end
mkClrDir(p.OutputDirectory);
%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...','Name',displFieldCorrProc.getName());
end

% Reading various constants
nFrames = movieData.nFrames_;

% Check displacement field process
iDisplFieldCalcProc =movieData.getProcessIndex('DisplacementFieldCalculationProcess',1,1);     
if isempty(iDisplFieldCalcProc)
    error(['Displacement field calculation has not been run! '...
        'Please run displacement field calculation prior to force field calculation!'])   
end

displFieldCalcProc=movieData.processes_{iDisplFieldCalcProc};
if ~displFieldCalcProc.checkChannelOutput
    error(['The channel must have a displacement field ! ' ...
        'Please calculate displacement field to all needed channels before '...
        'running force field calculation!'])
end
displParams = displFieldCalcProc.funParams_;
inFilePaths{1} = displFieldCalcProc.outFilePaths_{1};
displFieldCorrProc.setInFilePaths(inFilePaths);

% Set up the output directories
outputFile{1,1} = [p.OutputDirectory filesep 'displField.mat'];
outputFile{2,1} = [p.OutputDirectory filesep 'dispMaps.mat'];
mkClrDir(p.OutputDirectory);
displFieldCorrProc.setOutFilePaths(outputFile);

% get firstMask
iTFMPack = movieData.getPackageIndex('TFMPackage');
tfmPackageHere=movieData.packages_{iTFMPack}; iSDCProc=1;
SDCProc=tfmPackageHere.processes_{iSDCProc};
% iSDCProc =movieData.getProcessIndex('StageDriftCorrectionProcess',1,1);     
pDistProc = displFieldCalcProc.funParams_;
if ~isempty(SDCProc)
%     SDCProc=movieData.processes_{iSDCProc};
    if ~SDCProc.checkChannelOutput(pDistProc.ChannelIndex)
        error(['The channel must have been corrected ! ' ...
            'Please apply stage drift correction to all needed channels before '...
            'running displacement field calclation tracking!'])
    end
    refFrame = double(imread(SDCProc.outFilePaths_{2,pDistProc.ChannelIndex}));
else
    refFrame = double(imread(pDistProc.referenceFramePath));
end
firstMask=refFrame>0;
%% --------------- Displacement field correction ---------------%%% 

disp('Starting correcting displacement field...')
% Anonymous functions for reading input/output
displField=displFieldCalcProc.loadChannelOutput;


disp('Detecting and filtering vector field outliers...')
logMsg = 'Please wait, detecting and filtering vector field outliers';
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t/60)) 'min'];
tic;
useGrid=displParams.useGrid;

% %Parse input, store in parameter structure
% pd = parseProcessParams(displFieldCalcProc,paramsIn);

% Perform vector field outlier detection
if feature('ShowFigureWindows'), waitbar(0,wtBar,sprintf(logMsg)); end
if feature('ShowFigureWindows'),parfor_progress(nFrames); end

outlierThreshold=p.outlierThreshold;
% if useGrid
%     disp('In previous step, PIV was used, which does not require the current filtering step. skipping...')
% else
%     parfor j= 1:nFrames
    for j= 1:nFrames
        % Outlier detection
        dispMat = [displField(j).pos displField(j).vec];
        % Take out duplicate points (Sangyoon)
        [dispMat,~,~] = unique(dispMat,'rows'); %dispMat2 = dispMat(idata,:),dispMat = dispMat2(iudata,:)
        displField(j).pos=dispMat(:,1:2);
        displField(j).vec=dispMat(:,3:4);

        if ~isempty(outlierThreshold)
            [outlierIndex,sparselyLocatedIdx,~,neighborhood_distance(j)] = detectVectorFieldOutliersTFM(dispMat,outlierThreshold,1,'maxDist',2^(nextpow2(pDistProc.minCorLength)+1) + pDistProc.maxFlowSpeed);
            %displField(j).pos(outlierIndex,:)=[];
            %displField(j).vec(outlierIndex,:)=[];
            dispMat(outlierIndex,3:4)=NaN;
            dispMat(sparselyLocatedIdx,3:4)=NaN;
            % I deleted this part for later gap-closing
            % Filter out NaN from the initial data (but keep the index for the
            % outliers)
    %         ind= ~isnan(dispMat(:,3));
    %         dispMat=dispMat(ind,:);

            displField(j).pos=dispMat(:,1:2);
            displField(j).vec=dispMat(:,3:4);

            % I deleted this part because artificially interpolated vector can
            % cause more error or false force. - Sangyoon June 2013
    %         % Filling all NaNs with interpolated displacement vectors -
    %         % We also calculate the interpolated displacements with a bigger correlation length.
    %         % They are considered smoothed displacements at the data points. Sangyoon
    %         dispMat = [dispMat(:,2:-1:1) dispMat(:,2:-1:1)+dispMat(:,4:-1:3)];
    %         intDisp = vectorFieldSparseInterp(dispMat,...
    %             displField(j).pos(:,2:-1:1),...
    %             pd.minCorLength,pd.minCorLength,[],true);
    %         displField(j).vec = intDisp(:,4:-1:3) - intDisp(:,2:-1:1);
        end

        % Update the waitbar
    %     if mod(j,5)==1 && feature('ShowFigureWindows')
    %         tj=toc;
    %         waitbar(j/nFrames,wtBar,sprintf([logMsg timeMsg(tj*(nFrames-j)/j)]));
    %     end
        if feature('ShowFigureWindows'), parfor_progress; end
    end
% end
if feature('ShowFigureWindows'), parfor_progress(0); end

if p.fillVectors
    % Now this is the real cool step, to run trackStackFlow with known
    % information of existing displacement in neighbors
    % Check optional process Flow Tracking
    pStep2 = displParams;
    minCorLength = pStep2.minCorLength;
    enlargeFactor = 1;
    if ~isempty(SDCProc)
        s = load(SDCProc.outFilePaths_{3,pStep2.ChannelIndex},'T');
        residualT = s.T-round(s.T);
        refFrame = double(imread(SDCProc.outFilePaths_{2,pStep2.ChannelIndex}));
    else
        refFrame = double(imread(pStep2.referenceFramePath));
        residualT = zeros(nFrames,2);
    end
    logMsg = 'Please wait, retracking untracked points ...';
    timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t/60)) 'min'];
    tic

    % -- Pre-load all images to avoid repeated I/O inside parfor --
    allImages = cell(nFrames,1);
    for j = 1:nFrames
        if ~isempty(SDCProc)
            allImages{j} = double(SDCProc.loadChannelOutput(pStep2.ChannelIndex(1),j));
        else
            allImages{j} = double(movieData.channels_(pStep2.ChannelIndex(1)).loadImage(j));
        end
    end

    % -- Temporal warm-start: average displacement of adjacent frames --
    % For missing (p,j), use mean of tracked adjacent frames as initial estimate.
    % Dramatically reduces number of tracking iterations needed for time-lapse.
    nPointsAll = size(displField(1).pos,1);
    allVecXpre = nan(nPointsAll, nFrames);
    allVecYpre = nan(nPointsAll, nFrames);
    for j = 1:nFrames
        allVecXpre(:,j) = displField(j).vec(:,1);
        allVecYpre(:,j) = displField(j).vec(:,2);
    end
    temporalWarmWindow = min(3, max(1, floor(nFrames/2)));

    nFillingTries = 1000;
    unTrackedBeadsOrgAll = cell(nFrames,1);

    % -- WAVE PROPAGATION across frames (regular for loop) --
    % NOTE: parfor is intentionally NOT used here because
    % trackStackFlowWithHardCandidate already uses parfor internally.
    % Nested parfor would serialize the inner one, losing all parallelism.
    % The inner parfor (point-level) is more granular and gives better speedup.
    % Wave propagation: track most-constrained points first each iteration
    % -> geometric convergence vs random nFillingTries.
    for j = 1:nFrames
        dfj           = displField(j);
        currImage     = allImages{j};
        minNumNei     = 4;
        thresDist     = 4;
        nFailed       = 0;  nMaxFailed = 30;
        useNeighbors  = false;
        enlargeFactor = 1;
        minCorLengthJ = minCorLength;
        prevAttempt   = false;  prevNeiVecs = {};  prevIndices = [];

        unTrackedBeadsOrg = isnan(dfj.vec(:,1));
        unTrackedBeadsOrgAll{j} = unTrackedBeadsOrg;

        % Temporal warm-start: inject mean of adjacent frames for missing points
        tWin = max(1,j-temporalWarmWindow) : min(nFrames,j+temporalWarmWindow);
        tWin = tWin(tWin ~= j);
        if ~isempty(tWin)
            tempCandX = nanmean(allVecXpre(:, tWin), 2);
            tempCandY = nanmean(allVecYpre(:, tWin), 2);
            needsTempInit = unTrackedBeadsOrg & ~isnan(tempCandX);
            dfj.vec(needsTempInit,1) = tempCandX(needsTempInit);
            dfj.vec(needsTempInit,2) = tempCandY(needsTempInit);
        end

        nTracked = 1000;
        for k = 1:nFillingTries
            zeroDisp = dfj.vec(:,1)==0 & dfj.vec(:,2)==0;
            dfj.vec(zeroDisp,:) = NaN(sum(zeroDisp),2);
            unTrackedBeads = isnan(dfj.vec(:,1));
            ratioUntracked = sum(unTrackedBeads)/length(unTrackedBeads);

            if ratioUntracked < 0.0001 || (nTracked==0 && nFailed>nMaxFailed)
                break
            end

            currentBeads  = dfj.pos(unTrackedBeads,:);
            neighborBeads = dfj.pos(~unTrackedBeads,:);
            neighborVecs  = dfj.vec(~unTrackedBeads,:);

            % Single KDTree + binary search (FIX 1)
            maxSearchDist = neighborhood_distance(j);
            [idxFull, distFull] = KDTreeBallQuery(neighborBeads, currentBeads, maxSearchDist);
            lo = 0;  hi = maxSearchDist;
            for bsIter = 1:20
                mid = (lo + hi) / 2;
                if mean(cellfun(@(d) sum(d<=mid), distFull)) > minNumNei
                    hi = mid;
                else
                    lo = mid;
                end
            end
            optDist      = (lo + hi) / 2;
            idx          = cellfun(@(d) find(d<=optDist), distFull, 'UniformOutput', false);
            closeNeiVecs = cellfun(@(x) neighborVecs(x,:), idx, 'Unif', false);
            nNei         = cellfun(@numel, idx);
            idCloseEnough = nNei > thresDist+1;

            % O(1) prevLookup (FIX 3)
            if prevAttempt && ~isempty(prevIndices)
                nTotal = numel(unTrackedBeads);
                prevLookup = zeros(nTotal,1);
                prevLookup(prevIndices) = 1:numel(prevIndices);
                pp = 0;
                for jj = find(unTrackedBeads')
                    pp = pp + 1;
                    indInPrev = prevLookup(jj);
                    if indInPrev > 0
                        if abs(mean(closeNeiVecs{pp}(:,1)) - mean(prevNeiVecs{indInPrev}(:,1))) < 1e-2 && ...
                           abs(mean(closeNeiVecs{pp}(:,2)) - mean(prevNeiVecs{indInPrev}(:,2))) < 1e-2
                            idCloseEnough(pp) = false;
                        end
                    end
                end
            end

            v = NaN(sum(unTrackedBeads),2);
            if any(idCloseEnough)
                [v(idCloseEnough,:), nTracked] = trackStackFlowWithHardCandidate( ...
                    cat(3,refFrame,currImage), currentBeads(idCloseEnough,:), ...
                    minCorLengthJ, pStep2.minCorLength, 'maxSpd', pStep2.maxFlowSpeed, ...
                    'mode', pStep2.mode, 'hardCandidates', closeNeiVecs(idCloseEnough), ...
                    'magDiffThreshold', p.magDiffThreshold, ...
                    'angDiffThreshold', p.angDiffThreshold, ...
                    'enlargeFactor', enlargeFactor, 'useNeighbors', useNeighbors);
            else
                nTracked = 0;
            end

            dfj.vec(unTrackedBeads,:) = [v(:,1)+residualT(j,2)  v(:,2)+residualT(j,1)];
            prevAttempt = true;
            prevNeiVecs = closeNeiVecs;
            prevIndices = find(unTrackedBeads);

            if nTracked == 0
                nFailed      = nFailed + 1;
                useNeighbors = true;
                minNumNei    = 3 + nFailed;
                if minNumNei > 7 && minCorLengthJ > 7, minCorLengthJ = minCorLengthJ-2; end
                if enlargeFactor < 2 && nFailed > 1,   enlargeFactor = enlargeFactor*1.1; end
                if minNumNei > 7, minNumNei = 7; end
                enlargeFactor = 1;
            else
                nFailed = 0;
            end
        end % inner k loop

        displField(j) = dfj;
        disp(['Done for frame ' num2str(j) '/' num2str(nFrames) '.'])
        if feature('ShowFigureWindows')
            tj = toc;
            waitbar(j/nFrames, wtBar, sprintf([logMsg timeMsg(tj*(nFrames-j)/max(j,1))]));
        end
    end
    %Filtering again after fill-in
    if ~useGrid
        for j= 1:nFrames
            unTrackedBeadsOrg = unTrackedBeadsOrgAll{j};  % per-frame original mask
            dispMat = [displField(j).pos displField(j).vec];
            [dispMat,~,~] = unique(dispMat,'rows');
            displField(j).pos = dispMat(:,1:2);
            displField(j).vec = dispMat(:,3:4);
            [outlierIndex,sparselyLocatedIdx] = detectVectorFieldOutliersTFM(dispMat,outlierThreshold*10,1);
            outlierIndex       = outlierIndex(~ismember(outlierIndex',find(unTrackedBeadsOrg)')');
            sparselyLocatedIdx = sparselyLocatedIdx(~ismember(sparselyLocatedIdx',find(unTrackedBeadsOrg)')');
            dispMat(outlierIndex,3:4)       = NaN;
            dispMat(sparselyLocatedIdx,3:4) = NaN;
            displField(j).pos = dispMat(:,1:2);
            displField(j).vec = dispMat(:,3:4);
        end
    end
%     if feature('ShowFigureWindows'), parfor_progress(0); end
end
% Here, if nFrame>1, we do inter- and extrapolation of displacement vectors
% to prevent sudden, wrong force field change.
if nFrames>1 && ~displParams.useGrid
    disp('Performing displacement vector gap closing ...')
    % Depending on stage drift correction, some beads can be missed in certain
    % frames. Now it's time to make the same positions for all frames
    % go through each frame and filter points to the common ones in
    % iMinPointFrame - this needs to be improved by checking intersection
    % of all frames to find truly common beads, once there is error here.
    mostCommonPos = displField(1).pos;
    for ii= 2:nFrames
        commonPos=intersect(displField(ii).pos,mostCommonPos,'rows');
        mostCommonPos = commonPos;
    end
    for ii= 1:nFrames
        [commonPos,ia,~]=intersect(displField(ii).pos,mostCommonPos,'rows');
        displField(ii).pos = commonPos;
        displField(ii).vec = displField(ii).vec(ia,:);
    end
    % going through each point, see if there is NaN at each displacment
    % history and fill the gap
    logMsg = 'Performing displacement vector gap closing ...';

    nPoints = length(displField(1).pos(:,1));

    % FIX 2: Build full displacement matrices once (nPoints x nFrames).
    % Replaces per-point arrayfun(@(x) x.vec(k,1), displField) which
    % iterates over the whole struct array for every point k.
    allVecX = cell2mat(arrayfun(@(f) f.vec(:,1), displField(:)', 'UniformOutput', false)); % nPoints x nFrames
    allVecY = cell2mat(arrayfun(@(f) f.vec(:,2), displField(:)', 'UniformOutput', false));

    t = 1:nFrames;
    for k=1:nPoints
        curVecX = allVecX(k,:);
        curVecY = allVecY(k,:);
        nanMask = isnan(curVecX);
        if any(nanMask) && sum(~nanMask)/nFrames > 0.6
            t_nn = t(~nanMask);
            curVecX2 = interp1(t_nn, curVecX(~nanMask), t, 'linear');
            curVecY2 = interp1(t_nn, curVecY(~nanMask), t, 'linear');
            nanIdx = find(nanMask);
            for ii = nanIdx
                displField(ii).vec(k,:) = [curVecX2(ii) curVecY2(ii)];
            end
            allVecX(k, nanIdx) = curVecX2(nanIdx);
            allVecY(k, nanIdx) = curVecY2(nanIdx);
        end
        if mod(k,100)==1 && feature('ShowFigureWindows')
            tj=toc;
            waitbar(k/nPoints,wtBar,sprintf([logMsg timeMsg(tj*(nPoints-k)/k)]));
        end
    end
else
    if displParams.useGrid
        disp('In previous step, PIV was used, which does not require the current gap closing step. skipping...')
    end
end

% ================================================================
% Final filtering + 3D spatiotemporal filter (when nFrames > 1)
% ================================================================
% After fill-in, some re-tracked vectors may still be inaccurate.
% Run a final round of standard outlier detection first.
if p.fillVectors && ~useGrid
    disp('Final spatial outlier filtering after fill-in...')
    for j = 1:nFrames
        dispMat = [displField(j).pos displField(j).vec];
        [dispMat,~,~] = unique(dispMat,'rows');
        displField(j).pos = dispMat(:,1:2);
        displField(j).vec = dispMat(:,3:4);
        [outlierIndex, sparselyLocatedIdx] = detectVectorFieldOutliersTFM(dispMat, outlierThreshold, 1);
        % Preserve well-tracked original beads from being NaN'd
        outlierIndex       = outlierIndex(~ismember(outlierIndex', find(~unTrackedBeadsOrg)')');
        sparselyLocatedIdx = sparselyLocatedIdx(~ismember(sparselyLocatedIdx', find(~unTrackedBeadsOrg)')');
        dispMat(outlierIndex,    3:4) = NaN;
        dispMat(sparselyLocatedIdx, 3:4) = NaN;
        displField(j).pos = dispMat(:,1:2);
        displField(j).vec = dispMat(:,3:4);
    end
end

% 3D spatiotemporal filtering: uses spatial + temporal neighbors jointly
% for more robust outlier detection and gap-closing across the time-lapse.
if nFrames > 1 && ~displParams.useGrid
    disp('Running 3D spatiotemporal filter...')
    displField = filterDisplacementField3D(displField, ...
        'spatialRadius',    neighborhood_distance(1), ...
        'temporalWindow',   p.temporalWindow, ...
        'outlierThreshold', p.outlierThreshold3D, ...
        'spatialWeight',    p.spatialWeight3D, ...
        'temporalWeight',   p.temporalWeight3D, ...
        'verbose',          true);
end

% Find rotational registration
if p.doRotReg, displField=perfRotReg(displField); end 

%% Displacement map creation - this is shifted version
% Now decided to discard this part to save the hard disk space 20180421 SH
% [dMapIn, dmax, dmin, cropInfo,dMapXin,dMapYin,reg_grid] = generateHeatmapShifted(displField,displField,0);
% Insert displacement map in displField.pos 
% disp('Generating displacement maps ...')
% dMap = cell(1,nFrames);
% dMapX = cell(1,nFrames);
% dMapY = cell(1,nFrames);
% outputDir = fullfile(p.OutputDirectory,'displMaps');
% mkClrDir(outputDir);
% fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];
% numStr = @(frame) num2str(frame,fString);
% outFileDMap=@(frame) [outputDir filesep 'displMap' numStr(frame) '.mat'];
displFieldShifted(nFrames)=struct('pos','','vec','');
[reg_grid,~,~,~]=createRegGridFromDisplField(displField,1.0,0);

for ii=1:nFrames
    % starts with original size of beads
%     cur_dMap = zeros(size(firstMask));
%     cur_dMapX = zeros(size(firstMask));
%     cur_dMapY = zeros(size(firstMask));
%     cur_dMap(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = dMapIn{ii};
%     cur_dMapX(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = dMapXin{ii};
%     cur_dMapY(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = dMapYin{ii};
%     dMap{ii} = cur_dMap;
%     dMapX{ii} = cur_dMapX;
%     dMapY{ii} = cur_dMapY;
%     save(outFileDMap(ii),'cur_dMap','cur_dMapX','cur_dMapY'); % I removed v7.3 option to save the space,'-v7.3');
    % Shifted displField vector field
    [grid_mat,iu_mat, ~,~] = interp_vec2grid(displField(ii).pos, displField(ii).vec,[],reg_grid);
   
    [displFieldShiftedpos,displFieldShiftedvec, ~, ~] = interp_vec2grid(grid_mat+iu_mat, iu_mat,[],grid_mat); %1:cluster size
    pos = [reshape(displFieldShiftedpos(:,:,1),[],1) reshape(displFieldShiftedpos(:,:,2),[],1)]; %dense
    disp_vec = [reshape(displFieldShiftedvec(:,:,1),[],1) reshape(displFieldShiftedvec(:,:,2),[],1)]; 

    displFieldShifted(ii).pos = pos;
    displFieldShifted(ii).vec = disp_vec;
end
% clear dMapIn dMapXin dMapYin
disp('Saving ...')
save(outputFile{1},'displField','displFieldShifted','-v7.3');
% Saving the dMap which stores information
% dMap.eachDMapName = 'cur_dMap';
% dMap.outputDir = fullfile(p.OutputDirectory,'dislplMaps');
% dMap.outFileDMap = @(frame) [outputDir filesep 'displMap' numStr(frame) '.mat'];
dMap.displFieldPath = outputFile{1};
dMap.firstMaskSize = size(firstMask);

dMapX=dMap; %dMapX.eachTMapName = 'cur_dMapX';
dMapY=dMap; %dMapY.eachTMapName = 'cur_dMapY';
save(outputFile{2},'dMap','dMapX','dMapY'); % Updated, SH 20180225
% save(outputFile{2},'dMap','dMapX','dMapY','-v7.3'); % need to be updated for faster loading. SH 20141106
displMag=cell2mat(arrayfun(@(x) sqrt(x.vec(:,1).^2+x.vec(:,2).^2), displField,'unif',false)');
dmin = quantile(displMag(:),0.01); dmax = quantile(displMag(:),0.95);
displFieldCorrProc.setTractionMapLimits([dmin dmax])

%% Close waitbar
if feature('ShowFigureWindows'), close(wtBar); end

disp('Finished correcting displacement field!')