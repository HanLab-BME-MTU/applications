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
            [outlierIndex,sparselyLocatedIdx,~,neighborhood_distance(j)] = detectVectorFieldOutliersTFM(dispMat,outlierThreshold,1,'maxDist',pDistProc.minCorLength + pDistProc.maxFlowSpeed);
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
    nFillingTries=1000;
    for j= 1:nFrames
        % Read image and perform correlation
        if ~isempty(SDCProc)
            currImage = double(SDCProc.loadChannelOutput(pStep2.ChannelIndex(1),j));
        else
            currImage = double(movieData.channels_(pStep2.ChannelIndex(1)).loadImage(j));
        end
        k2=0;
        nTracked=1000; 
        useNeighbors = false;
        nFailed=0; nMaxFailed=30;
        prevAttempt=false; prevNeiVecs=[]; prevIndices=[];
        minNumNei=4; %was 3 before. This means there should be at least thresDist neighboring vectors per position. 
        % This is used for determining search radius. This will increase
        % over iteration.
        thresDist = 4; %this is used for determining what points will be interrogated.
        iFigureDrawing=0;
        unTrackedBeadsOrg = isnan(displField(j).vec(:,1));
        
        for k=1:nFillingTries
            iFigureDrawing=iFigureDrawing+1;
            % only un-tracked vectors, and ones with zero displacement (
            % I don't know where the zero displacements came from, but we
            % need to make these to NaNs)
            %Converting zero vectors into NaNs (It just make no sence that the
            %absolute zero displacement)
            zeroDisp = displField(j).vec(:,1)==0 & displField(j).vec(:,2)==0;
            displField(j).vec(zeroDisp,:) = NaN(sum(zeroDisp),2);
            unTrackedBeads=isnan(displField(j).vec(:,1));
            ratioUntracked = sum(unTrackedBeads)/length(unTrackedBeads);
            if ratioUntracked<0.0001 || (nTracked==0 && nFailed>nMaxFailed)
                break
            end
            currentBeads = displField(j).pos(unTrackedBeads,:);
            neighborBeads = displField(j).pos(~unTrackedBeads,:);
            neighborVecs = displField(j).vec(~unTrackedBeads,:);
            % Get neighboring vectors from these vectors (meanNeiVecs)
%             [idx] = KDTreeBallQuery(neighborBeads, currentBeads, (1-5*k/nFillingTries)*neighborhood_distance(j)); % Increasing search radius with further iteration
%             [idx,dist] = KDTreeClosestPoint(neighborBeads, currentBeads); % Increasing search radius with further iteration
            [idx] = KDTreeBallQuery(neighborBeads, currentBeads, (1)*neighborhood_distance(j)); % Decreasing search radius with further iteration
            % Here we want to be very wise. The fact that one bead location
            % has few  neighboring vectors is not always good. What if the
            % neighboring vectors are far from the location? 
            % We need to get an optimal distance 
            % First get the number distribution
            quantityNeis = cellfun(@numel,idx);
            qq=1;
            meanNumNeis = NaN(51,1);
            meanNumNeis(qq,1) = mean(quantityNeis);
            % Currently mean is about 30. I am going now keep lowering the
            % distance and see how mean is changing
            for k2=0.1:0.1:5
                [idx] = KDTreeBallQuery(neighborBeads, currentBeads, (1-0.2*k2)*neighborhood_distance(j)); % Decreasing search radius with further iteration
                quantityNeis = cellfun(@numel,idx);
                qq=qq+1;
                meanNumNeis(qq,1) = mean(quantityNeis);
            end
            %figure, plot(0:0.1:5',meanNumNeis)
            % Now limit the distance so that meanNumNei becomes only
            % thresDist
            [~,indOptDist]=min(abs(meanNumNeis -minNumNei));
            optDist = (1-0.015*indOptDist)*neighborhood_distance(j);
            [idx] = KDTreeBallQuery(neighborBeads, currentBeads, optDist); % optimal neighbors
            
%             % In case of empty idx, search with larger radius.
%             emptyCases = cellfun(@isempty,idx);
%             mulFactor=1;
%             while any(emptyCases)
%                 mulFactor=mulFactor+0.5;
%                 idxEmpty = KDTreeBallQuery(neighborBeads, currentBeads(emptyCases,:), mulFactor*(1+5*k/nFillingTries)*neighborhood_distance(j));
%                 idx(emptyCases)=idxEmpty;
%                 emptyCases = cellfun(@isempty,idx);
%             end
            % Subsample idx to reduce computing time
%             % Calculate the subsampling rate
%             leap = cellfun(@(x) max(1,round(length(x)/100)),idx,'Unif',false);
%             idx = cellfun(@(x,y) x(1:y:end,1),idx,leap,'Unif',false);
            closeNeiVecs = cellfun(@(x) neighborVecs(x,:),idx,'Unif',false);
%             closeNeiVecs = arrayfun(@(x) neighborVecs(x,:),idx,'Unif',false);
%             idCloseEnough = dist<thresDist;
            atLeast5neis = cellfun(@(x) numel(x)>thresDist+1,idx);
            idCloseEnough = atLeast5neis; %dist<thresDist;
            % See if there is any previous attempt
            tolR = 1e-2;
            if prevAttempt
                % Collect bead locations that have no neighbor information
                % change.
                curUntrackedIndeces = 1;
                pp=0;
                for jj=find(unTrackedBeads')
                    % We compare with mean vec
                    pp=pp+1;
                    if ismember(jj,prevIndices)
                        % Find the index inside prevNeiVecs
                        indInPrevNeiVecs = find(jj==prevIndices);
                        curMeanX = mean(closeNeiVecs{pp}(:,1));
                        curMeanY = mean(closeNeiVecs{pp}(:,2));
                        prevMeanX = mean(prevNeiVecs{indInPrevNeiVecs}(:,1));
                        prevMeanY = mean(prevNeiVecs{indInPrevNeiVecs}(:,2));
                        if abs(curMeanX-prevMeanX)<tolR && abs(curMeanY-prevMeanY)<tolR
                            idCloseEnough(pp)=false;
                        end
                    end
                end
            end
            
            v=NaN(sum(unTrackedBeads),2);
        %     meanNeiVecs = cellfun(@mean,closeNeiVecs,'Unif',false);

            [v(idCloseEnough,:),nTracked] = trackStackFlowWithHardCandidate(cat(3,refFrame,currImage),currentBeads(idCloseEnough,:),...
                minCorLength,pStep2.minCorLength,'maxSpd',pStep2.maxFlowSpeed,...
                'mode',pStep2.mode,'hardCandidates',closeNeiVecs(idCloseEnough),'magDiffThreshold',p.magDiffThreshold,...
                'angDiffThreshold',p.angDiffThreshold,'enlargeFactor',enlargeFactor,...
                'useNeighbors',useNeighbors); %, 'hardCandidateDists', dist(idCloseEnough));%,'usePIVSuite', pStep2.usePIVSuite);

%             if j== 1 && (iFigureDrawing==1 || nFailed==nMaxFailed)
            if k==1 % This is for publication purpose for PTVR
                figure
            end
            if nansum(v(idCloseEnough,1))>0
                hold off
                quiver(neighborBeads(:,1),neighborBeads(:,2),neighborVecs(:,1),neighborVecs(:,2),0,'k')
                hold on
                curTrackedPos = ~isnan(v(:,1));
                plot(currentBeads(~curTrackedPos,1),currentBeads(~curTrackedPos,2),'bo')
                quiver(currentBeads(curTrackedPos,1),currentBeads(curTrackedPos,2),...
                    v(curTrackedPos,1)+residualT(j,2), v(curTrackedPos,2)+residualT(j,1),0,'r')
                set(gca, 'YDir','reverse')
                set(gca, 'XLim',[32 480], 'YLim',[32 480])
                set(gca, 'PlotBoxAspectRatio',[1 1 1])
                set(gca, 'XTickLabel', [])
                set(gca, 'YTickLabel', [])
                savefig([p.OutputDirectory filesep 'field_iter' num2str(k) '.fig'])
                print([p.OutputDirectory filesep 'field_iter' num2str(k) '.tif'], '-djpeg')
            end
%             end
        %     displField(j).pos(unTrackedBeads,:)=currentBeads; % validV is removed to include NaN location - SH 030417
            displField(j).vec(unTrackedBeads,:)=[v(:,1)+residualT(j,2) v(:,2)+residualT(j,1)]; % residual should be added with oppiste order! -SH 072514

%             % Filtering it again because retracked vectors might contain
%             % errors - I decided to not do this but controlling the
%             significance criteria
%             dispMat = [displField(j).pos displField(j).vec];
%             [outlierIndex,sparselyLocatedIdx,~,neighborhood_distance(j)] = detectVectorFieldOutliersTFM(dispMat,outlierThreshold,1);
%             if iFigureDrawing==1 || nTracked==nMaxFailed
%                 quiver(dispMat(outlierIndex,1),dispMat(outlierIndex,2),dispMat(outlierIndex,3), dispMat(outlierIndex,4),0,'g')
%                 quiver(dispMat(sparselyLocatedIdx,1),dispMat(sparselyLocatedIdx,2),dispMat(sparselyLocatedIdx,3), dispMat(sparselyLocatedIdx,4),0,'b')
%             end
%             dispMat(outlierIndex,3:4)=NaN;
%             dispMat(sparselyLocatedIdx,3:4)=NaN;
%             displField(j).pos=dispMat(:,1:2);
%             displField(j).vec=dispMat(:,3:4);
%             nTracked = sum(~isnan(dispMat(:,3)) & unTrackedBeads);
%             disp(['Tracked - filtered = ' num2str(nTracked)])

            % We don't need to track ones that have failed tracking if
            % the radius is not updated (thus the neighboring
            % information is not updated). 
            % Save untracked locations
            prevAttempt=true;
            prevNeiVecs = closeNeiVecs;
            prevIndices = find(unTrackedBeads);

            if nTracked==0
                nFailed=nFailed+1;
                useNeighbors = true;
                minNumNei = 3 +nFailed; %(-1)^nFailed*ceil((nFailed+1)/2);
                if minNumNei > 7 && minCorLength > 7
                    minCorLength = minCorLength - 2;
                end
                if enlargeFactor<2 && nFailed>1 % perform this after all median-based vectors are found
                    enlargeFactor = enlargeFactor *1.1;
                end
%                 if thresDist==0
%                     thresDist=nFailed+1;
%                 end
            else
                nFailed=0;
                useNeighbors = false;
                if minCorLength < 8
                    minCorLength = pStep2.minCorLength;
                end
                if minNumNei > 7
                    minNumNei = 7; % Decided to let thresDist stay where it was, but there is a maximum.
                end
                enlargeFactor=1;
            end
        end
        disp(['Done for frame ' num2str(j) '/' num2str(nFrames) '.'])
        % Update the waitbar
        if feature('ShowFigureWindows')
            tj=toc;
            waitbar(j/nFrames,wtBar,sprintf([logMsg timeMsg(tj*(nFrames-j)/j)]));
        end
    end
    %Filtering again
    if ~useGrid
        for j= 1:nFrames
            % Outlier detection
            dispMat = [displField(j).pos displField(j).vec];
            % Take out duplicate points (Sangyoon)
            [dispMat,~,~] = unique(dispMat,'rows'); %dispMat2 = dispMat(idata,:),dispMat = dispMat2(iudata,:)
            displField(j).pos=dispMat(:,1:2);
            displField(j).vec=dispMat(:,3:4);

            [outlierIndex,sparselyLocatedIdx] = detectVectorFieldOutliersTFM(dispMat,outlierThreshold*10,1);
            %displField(j).pos(outlierIndex,:)=[];
            %displField(j).vec(outlierIndex,:)=[];
            % We neglect these if the outliers are from well-tracked ones:
            outlierIndex = outlierIndex(~ismember(outlierIndex',find(unTrackedBeadsOrg)')');
            sparselyLocatedIdx = sparselyLocatedIdx(~ismember(sparselyLocatedIdx',find(unTrackedBeadsOrg)')');
            dispMat(outlierIndex,3:4)=NaN;
            dispMat(sparselyLocatedIdx,3:4)=NaN;

            displField(j).pos=dispMat(:,1:2);
            displField(j).vec=dispMat(:,3:4);
%             if feature('ShowFigureWindows'), parfor_progress; end
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
    for k=1:nPoints
        % build each disp vector history
        curVecX = arrayfun(@(x) x.vec(k,1),displField);
        curVecY = arrayfun(@(x) x.vec(k,2),displField);
        if any(isnan(curVecX)) && sum(~isnan(curVecX))/nFrames>0.6
            t = 1:length(curVecX);
            t_nn = t(~isnan(curVecX));
            curVecX2 = interp1(t_nn,curVecX(~isnan(curVecX)),t,'linear');
            curVecY2 = interp1(t_nn,curVecY(~isnan(curVecX)),t,'linear');
            for ii=find(isnan(curVecX))
                displField(ii).vec(k,:) = [curVecX2(ii) curVecY2(ii)];
            end
        else
            continue
        end
        if mod(k,5)==1 && feature('ShowFigureWindows')
            tj=toc;
            waitbar(k/nPoints,wtBar,sprintf([logMsg timeMsg(tj*(nPoints-k)/k)]));
        end
    end
else
    if displParams.useGrid
        disp('In previous step, PIV was used, which does not require the current gap closing step. skipping...')
    end
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