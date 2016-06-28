function calculateMovieDisplacementField(movieData,varargin)
% calculateMovieDisplacementField calculate the displacement field
%
% calculateMovieDisplacementField 
%
% SYNOPSIS calculateMovieDisplacementField(movieData,paramsIn)
%
% INPUT   
%   movieData - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%

% Sebastien Besson, Sep 2011

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous stage drift processes                                                                     
iProc = movieData.getProcessIndex('DisplacementFieldCalculationProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(DisplacementFieldCalculationProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end
displFieldProc = movieData.processes_{iProc};
%Parse input, store in parameter structure
p = parseProcessParams(displFieldProc,paramsIn);
addNonLocMaxBeads = p.addNonLocMaxBeads;
%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...','Name',displFieldProc.getName());
else
    wtBar = -1;
end

% Reading various constants
nFrames = movieData.nFrames_;

% Check optional process Flow Tracking
iSDCProc =movieData.getProcessIndex('StageDriftCorrectionProcess',1,1);     
if ~isempty(iSDCProc)
    SDCProc=movieData.processes_{iSDCProc};
    if ~SDCProc.checkChannelOutput(p.ChannelIndex)
        error(['The channel must have been corrected ! ' ...
            'Please apply stage drift correction to all needed channels before '...
            'running displacement field calclation tracking!'])
    end
    imDirs{1} = SDCProc.outFilePaths_{1,p.ChannelIndex};
    s = load(SDCProc.outFilePaths_{3,p.ChannelIndex},'T');
    residualT = s.T-round(s.T);
    refFrame = double(imread(SDCProc.outFilePaths_{2,p.ChannelIndex}));
else
    imDirs  = movieData.getChannelPaths(p.ChannelIndex);
    refFrame = double(imread(p.referenceFramePath));
    residualT = zeros(nFrames,2);
end
inFilePaths{1,p.ChannelIndex} = imDirs{:};
displFieldProc.setInFilePaths(inFilePaths);
    
% Set up the output directories
outputFile{1,1} = [p.OutputDirectory filesep 'displField.mat'];
outputFile{2,1} = [p.OutputDirectory filesep 'dispMaps.mat'];

% Add a recovery mechanism if process has been stopped in the middle of the
% computation to re-use previous results
firstFrame =1; % Set the strating fram eto 1 by default
if exist(outputFile{1},'file');
    % Check analyzed frames
    s=load(outputFile{1},'displField');
    frameDisplField=~arrayfun(@(x)isempty(x.pos),s.displField);
    
    if ~all(frameDisplField) && ~all(~frameDisplField) && usejava('desktop')
        % Look at the first non-analyzed frame
        firstFrame = find(~frameDisplField,1);
        % Ask the user if display mode is active
        if ishandle(wtBar),
            recoverRun = questdlg(...
                ['A displacement field output has been dectected with ' ...
                num2str(firstFrame-1) ' analyzed frames. Do you' ...
                ' want to use these results and continue the analysis'],...
                'Recover previous run','Yes','No','Yes');
            if ~strcmpi(recoverRun,'Yes'), firstFrame=1; end
        end
    end
end

if firstFrame == 1, 
    % Backup the original vectors to backup folder
    display('Backing up the original data')
    backupFolder = [p.OutputDirectory ' Backup']; % name]);
    if exist(p.OutputDirectory,'dir')
        ii = 1;
        while exist(backupFolder,'dir')
            backupFolder = [p.OutputDirectory ' Backup ' num2str(ii)];
            ii=ii+1;
        end
        mkdir(backupFolder);
        copyfile(p.OutputDirectory, backupFolder)
    end

    % Clean output file and initialize displacement field structure
    mkClrDir(p.OutputDirectory); 
    displField(nFrames)=struct('pos',[],'vec',[]);
else
    % Load old displacement field structure 
    displField=s.displField;
end

displFieldProc.setOutFilePaths(outputFile);

%% --------------- Displacement field calculation ---------------%%% 

disp('Starting calculating displacement field...')
% Get the mask
maskArray = movieData.getROIMask;
% Use mask of first frame to filter bead detection
firstMask = refFrame>0; %false(size(refFrame));
tempMask = maskArray(:,:,1);
% firstMask(1:size(tempMask,1),1:size(tempMask,2)) = tempMask;
tempMask2 = false(size(refFrame));
y_shift = find(any(firstMask,2),1);
x_shift = find(any(firstMask,1),1);

tempMask2(y_shift:y_shift+size(tempMask,1)-1,x_shift:x_shift+size(tempMask,2)-1) = tempMask;
firstMask = tempMask2 & firstMask;

    % if ~p.useGrid
% end

% if p.useGrid && p.highRes
%     tempDisplField.pos = beads;
%     [~,xvec,yvec,~]=createRegGridFromDisplField(tempDisplField,2,0);
%     beads = [xvec yvec];
% end

disp('Calculating displacement field...')
logMsg = 'Please wait, calculating displacement field';
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t/60)) 'min'];
tic;

% Perform sub-pixel registration
if ishandle(wtBar), waitbar(0,wtBar,sprintf(logMsg)); end

for j= firstFrame:nFrames
    % Read image and perform correlation
    if ~isempty(iSDCProc)
        currImage = double(SDCProc.loadChannelOutput(p.ChannelIndex(1),j));
    else
        currImage = double(movieData.channels_(p.ChannelIndex(1)).loadImage(j));
    end
    if ~p.useGrid
    % if strcmp(movieData.getChannel(p.ChannelIndex).imageType_,'Widefield')
        % Detect beads in reference frame
        disp('Determining PSF sigma from reference frame...')
        % Adaptation of psfSigma from bead channel image data
        psfSigma = getGaussianSmallestPSFsigmaFromData(refFrame,'Display',false);
        if isnan(psfSigma) || psfSigma>movieData.channels_(1).psfSigma_*3 
            if strcmp(movieData.getChannel(p.ChannelIndex(1)).imageType_,'Widefield') || movieData.pixelSize_>130
                psfSigma = movieData.channels_(1).psfSigma_*2; %*2 scale up for widefield
            elseif strcmp(movieData.getChannel(p.ChannelIndex(1)).imageType_,'Confocal')
                psfSigma = movieData.channels_(1).psfSigma_*0.79; %*4/7 scale down for  Confocal finer detection SH012913
            elseif strcmp(movieData.getChannel(p.ChannelIndex(1)).imageType_,'TIRF')
                psfSigma = movieData.channels_(1).psfSigma_*3/7; %*3/7 scale down for TIRF finer detection SH012913
            else
                error('image type should be chosen among Widefield, confocal and TIRF!');
            end
        end
        disp(['Determined sigma: ' num2str(psfSigma)])

        disp('Detecting beads in the reference frame...')
        pstruct = pointSourceDetection(refFrame, psfSigma, 'alpha', p.alpha,'Mask',firstMask,'FitMixtures',true);
        assert(~isempty(pstruct), 'Could not detect any bead in the reference frame');
        % filtering out points in saturated image based on pstruct.c
        [N,edges]= histcounts(pstruct.c);
        % starting with median, find a edge disconnected with two consequtive
        % zeros.
        medC = median(pstruct.c);
        idxAfterMedC=find(edges>medC);
        qq=idxAfterMedC(1);
        while N(qq)>0 || N(qq+1)>0
            qq=qq+1;
            if qq>=length(edges)-1
                break
            end
        end
        idx = pstruct.c<edges(qq);
        beads = [round(pstruct.x(idx)') round(pstruct.y(idx)')];
        %     beads = [ceil(pstruct.x') ceil(pstruct.y')];

        % Subsample detected beads ensuring beads are separated by at least half of
        % the correlation length - commented out to get more beads
        if ~p.highRes
            disp('Subsampling detected beads (normal resolution)...')
            max_beads_distance = floor(p.minCorLength/2);
        else
            % To get high-resolution information, subsample detected beads ensuring 
            % beads are separated by 0.1 um the correlation length 
            disp('Subsampling detected beads (high resolution)...')
            max_beads_distance = (100/movieData.pixelSize_);
        end

        idx = KDTreeBallQuery(beads, beads, max_beads_distance);
        valid = true(numel(idx),1);
        for i = 1 : numel(idx)
            if ~valid(i), continue; end
            neighbors = idx{i}(idx{i} ~= i);
            valid(neighbors) = false;
        end
        beads = beads(valid, :);
        % It doesn't critically require local maximal pixel to start
        % x-correlation-based tracking. Thus, to increase spatial resolution,
        % we add additional points in the mid points of pstruct
        % We first randomly distribute point, and if it is not too close to
        % existing points and the intensity of the point is above a half of the
        % existing points, include the point into the point set
        if addNonLocMaxBeads
            disp('Finding additional non-local-maximal points with high intensity ...')
            distance=zeros(length(beads),1);
            for i=1:length(beads)
                neiBeads = beads;
                neiBeads(i,:)=[];
                [~,distance(i)] = KDTreeClosestPoint(neiBeads,beads(i,:));
            end
            avg_beads_distance = quantile(distance,0.4);%mean(distance);%size(refFrame,1)*size(refFrame,2)/length(beads);
            notSaturated = true;
            xmin = min(pstruct.x);
            xmax = max(pstruct.x);
            ymin = min(pstruct.y);
            ymax = max(pstruct.y);
        %     avgAmp = mean(pstruct.A);
        %     avgBgd = mean(pstruct.c);
        %     thresInten = avgBgd+0.02*avgAmp;
            thresInten = quantile(pstruct.c,0.1);
            maxNumNotDetected = 100; % the number of maximum trial without detecting possible point
            numNotDetected = 0;
            numPrevBeads = size(beads,1);
            while notSaturated
                x_new = xmin + (xmax-xmin)*rand();
                y_new = ymin + (ymax-ymin)*rand();
                [~,distToPoints] = KDTreeClosestPoint(beads,[x_new,y_new]);
                inten_new = refFrame(round(y_new),round(x_new));
                if distToPoints>avg_beads_distance && inten_new>thresInten
                    beads = [beads; x_new, y_new];
                    numNotDetected = 0;
                else
                    numNotDetected=numNotDetected+1;
                end
                if numNotDetected>maxNumNotDetected
                    notSaturated = false; % this means now we have all points to start tracking from the image
                end
            end
            disp([num2str(size(beads,1)-numPrevBeads) ' points were additionally detected for fine tracking. Total detected beads: ' num2str(length(beads))])
        end
        % Exclude all beads which are less  than half the correlation length 
        % away from the padded border. By default, no centered template should 
        % include any NaN's for correlation
        % Create beads mask with zero intensity points as false
        beadsMask = true(size(refFrame));
        % beadsMask(currImage==0)=false;
        % Remove false regions non-adjacent to the image border
        beadsMask = beadsMask | imclearborder(~beadsMask);
        %         % Erode the mask with half the correlation length and filter beads
        %         erosionDist=round((p.minCorLength+1)/2);
        % Erode the mask with the correlation length + half maxFlowSpeed
        % and filter beads to minimize error
        if p.noFlowOutwardOnBorder
            erosionDist=(p.minCorLength+1);
        else
            erosionDist=p.minCorLength+1+round(p.maxFlowSpeed/4);
        end            
        beadsMask=bwmorph(beadsMask,'erode',erosionDist);
        %         beadsMask=imerode(beadsMask,strel('square',erosionDist));
        indx=beadsMask(sub2ind(size(beadsMask),ceil(beads(:,2)), ceil(beads(:,1))));
        localbeads = beads(indx,:);

        %     % Select only beads which are min correlation length away from the border of the
        %     % reference frame 
        %     beadsMask = true(size(refFrame));
        %     erosionDist=p.minCorLength;
        %     % erosionDist=p.minCorLength+1;
        %     beadsMask(erosionDist:end-erosionDist,erosionDist:end-erosionDist)=false;
        %     indx=beadsMask(sub2ind(size(beadsMask),ceil(beads(:,2)),ceil(beads(:,1))));
        %     beads(indx,:)=[];

        % Track beads displacement in the xy coordinate system
        v = trackStackFlow(cat(3,refFrame,currImage),localbeads,...
            p.minCorLength,p.minCorLength,'maxSpd',p.maxFlowSpeed,...
            'mode',p.mode);

        % Extract finite displacement and prepare displField structure in the xy
        % coordinate system
        validV = ~isinf(v(:,1));
        displField(j).pos=localbeads(validV,:);
        displField(j).vec=[v(validV,1)+residualT(j,2) v(validV,2)+residualT(j,1)]; % residual should be added with oppiste order! -SH 072514
    else
        pivPar = [];      % variable for settings
        pivData = [];     % variable for storing results

        [pivPar, pivData] = pivParams(pivData,pivPar,'defaults');     
        % Set the size of interrogation areas via fields |iaSizeX| and |iaSizeY| of |pivPar| variable:
%         pivPar.iaSizeX = [64 32 16 2^(nextpow2(p.minCorLength)-1)];     % size of interrogation area in X 
        nextPow2=nextpow2(p.minCorLength);
        BiggestSize=2^(nextPow2+1);
        SecondSize=2^(nextPow2);
        ThirdSize=2^(nextPow2-1);
        FourthSize=2^(nextPow2-2);
        pivPar.iaSizeX = [BiggestSize SecondSize ThirdSize ThirdSize];     % size of interrogation area in X 
        pivPar.iaStepX = [SecondSize SecondSize ThirdSize FourthSize];     % grid spacing of velocity vectors in X
        pivPar.iaSizeY = [BiggestSize SecondSize ThirdSize ThirdSize];     % size of interrogation area in X 
        pivPar.iaStepY = [SecondSize SecondSize ThirdSize FourthSize];    % grid spacing of velocity vectors in X
%         pivPar.iaSizeX = [64 32 16 8];     % size of interrogation area in X 
%         pivPar.iaStepX = [32 16  8 4];     % grid spacing of velocity vectors in X
%         pivPar.iaSizeY = [64 32 16 8];     % size of interrogation area in X 
%         pivPar.iaStepY = [32 16  8 4];     % grid spacing of velocity vectors in X
% %         pivPar.iaStepX = [32 16  8 8];     % grid spacing of velocity vectors in X
% %         pivPar.iaSizeY = [64 32 16 16];     % size of interrogation area in Y 
% %         pivPar.iaStepY = [32 16  8 8];     % grid spacing of velocity vectors in Y
        pivPar.ccWindow = 'Gauss2';   % This filter is relatively narrow and will 
        pivPar.smMethod = 'none';
        pivPar.iaMethod = 'defspline';
        pivPar.iaImageInterpolationMethod = 'spline';
        pivPar.imMask1=firstMask;
        pivPar.imMask2=firstMask;
%         pivData.X=localbeads(:,1);
%         pivData.Y=localbeads(:,2);
%         pivData.U=zeros(size(pivData.X));
%         pivData.V=zeros(size(pivData.Y));

        [pivData] = pivAnalyzeImagePair(refFrame,currImage,pivData,pivPar);
        validV = ~isnan(pivData.V);
        
        displField(j).pos=[pivData.X(validV), pivData.Y(validV)];
        displField(j).vec=[pivData.U(validV)+residualT(j,2), pivData.V(validV)+residualT(j,1)]; % residual should be added with oppiste order! -SH 072514

%         % testing additional pass of piv processing
%         pivPar.iaSizeX=[ 8];
%         pivPar.iaStepX=[ 4];
%         pivPar.anVelocityEst = 'pivData'; % use velocity data stored in pivData as velocity estimate used for image deformation. 
%         % By this setting, results of previous passes are transferred
%         pivPar.ccMethod = 'dcn';
%         pivPar.qvPair = {'U', 'clipLo', -0.2, 'clipHi', 0.2};
%         figure;
%         [pivPar2] = pivParams([],pivPar,'defaults');
%         pivData2 = pivAnalyzeImagePair(refFrame,currImage,pivData,pivPar2);
%         
%         displField(j).pos=[pivData2.X(:), pivData2.Y(:)];
%         displField(j).vec=[pivData2.U(:)+residualT(j,2), pivData2.V(:)+residualT(j,1)]; % residual should be added with oppiste order! -SH 072514
    end
    
    % Update the waitbar
    if mod(j,5)==1 && ishandle(wtBar)
        tj=toc;
        waitbar(j/nFrames,wtBar,sprintf([logMsg ...
            timeMsg(tj*(nFrames-firstFrame++1-j)/j)]));
    end
    
    % Save each iteration (for recovery of unfinished processes)
    save(outputFile{1},'displField');
end
%% Displacement map creation - this is shifted version
[dMapIn, dmax, dmin, cropInfo,dMapXin,dMapYin,reg_grid] = generateHeatmapShifted(displField,displField,0);
display(['Estimated displacement maximum = ' num2str(dmax) ' pixel.'])
% Insert traction map in forceField.pos 
disp('Generating displacement maps ...')
dMap = cell(1,nFrames);
dMapX = cell(1,nFrames);
dMapY = cell(1,nFrames);
displFieldShifted(nFrames)=struct('pos','','vec','');
for ii=1:nFrames
    % starts with original size of beads
    cur_dMap = zeros(size(firstMask));
    cur_dMapX = zeros(size(firstMask));
    cur_dMapY = zeros(size(firstMask));
    cur_dMap(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = dMapIn{ii};
    cur_dMapX(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = dMapXin{ii};
    cur_dMapY(cropInfo(2):cropInfo(4),cropInfo(1):cropInfo(3)) = dMapYin{ii};
    dMap{ii} = cur_dMap;
    dMapX{ii} = cur_dMapX;
    dMapY{ii} = cur_dMapY;
    % Shifted displField vector field
    [grid_mat,iu_mat, ~,~] = interp_vec2grid(displField(ii).pos, displField(ii).vec,[],reg_grid);
   
    [displFieldShiftedpos,displFieldShiftedvec, ~, ~] = interp_vec2grid(grid_mat+iu_mat, iu_mat,[],grid_mat); %1:cluster size
    pos = [reshape(displFieldShiftedpos(:,:,1),[],1) reshape(displFieldShiftedpos(:,:,2),[],1)]; %dense
    disp_vec = [reshape(displFieldShiftedvec(:,:,1),[],1) reshape(displFieldShiftedvec(:,:,2),[],1)]; 

    displFieldShifted(ii).pos = pos;
    displFieldShifted(ii).vec = disp_vec;
end
disp('Saving ...')
save(outputFile{1},'displField','displFieldShifted','-v7.3');
save(outputFile{2},'dMap','dMapX','dMapY','-v7.3'); % need to be updated for faster loading. SH 20141106
displFieldProc.setTractionMapLimits([dmin dmax])

% Close waitbar
if ishandle(wtBar), close(wtBar); end

disp('Finished calculating displacement field!')