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
% Sangyoon Han, From Oct 2012. Last updated: May 2017

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
p.usePIVSuite=true;
%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...','Name',displFieldProc.getName());
else
    wtBar = -1;
end

% Reading various constants
nFrames = movieData.nFrames_;

% Check optional process Flow Tracking
% iSDCProc = movieData.getProcessIndex('StageDriftCorrectionProcess',1,1);     
iTFMPack = movieData.getPackageIndex('TFMPackage');
tfmPackageHere=movieData.packages_{iTFMPack}; iSDCProc=1;
SDCProc=tfmPackageHere.processes_{iSDCProc};
if ~isempty(SDCProc)
%     SDCProc=movieData.processes_{iSDCProc};
    if ~SDCProc.checkChannelOutput(p.ChannelIndex)
        error(['The channel must have been corrected ! ' ...
            'Please apply stage drift correction to all needed channels before '...
            'running displacement field calclation tracking!'])
    end
    imDirs{1} = SDCProc.outFilePaths_{1,p.ChannelIndex};
    s = load(SDCProc.outFilePaths_{3,p.ChannelIndex},'T');
    residualT = s.T-round(s.T);
    T = s.T;
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
if exist(outputFile{1},'file')
    % Check analyzed frames
    sDisp=load(outputFile{1},'displField');
    frameDisplField=~arrayfun(@(x)isempty(x.pos),sDisp.displField);
    
    if ~all(frameDisplField) && ~all(~frameDisplField) && usejava('desktop')
        % Look at the first non-analyzed frame
        firstFrame = find(~frameDisplField,1);
        % Ask the user if display mode is active
        if ishandle(wtBar)
            recoverRun = questdlg(...
                ['A displacement field output has been dectected with ' ...
                num2str(firstFrame-1) ' analyzed frames. Do you' ...
                ' want to use these results and continue the analysis'],...
                'Recover previous run','Yes','No','Yes');
            if ~strcmpi(recoverRun,'Yes'), firstFrame=1; end
        end
    end
end

if firstFrame == 1 
    % Backup the original vectors to backup folder
    disp('Backing up the original data')
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
    displField=sDisp.displField;
end

displFieldProc.setOutFilePaths(outputFile);

%% --------------- Displacement field calculation ---------------%%% 

disp('Starting calculating displacement field...')
% Get the mask
maskArray = movieData.getROIMask;
% Use mask of first frame to filter bead detection
if ~isempty(SDCProc)
    firstImage = double(SDCProc.loadChannelOutput(p.ChannelIndex(1),1));
else
    firstImage = double(movieData.channels_(p.ChannelIndex(1)).loadImage(1));
end
firstMask = refFrame>0 & firstImage>0; %false(size(refFrame));
tempMask = maskArray(:,:,1);
% firstMask(1:size(tempMask,1),1:size(tempMask,2)) = tempMask;
ref_obj = imref2d(size(refFrame));
% tempMask2 = false(size(refFrame));
if ~isempty(SDCProc)
    if isa(SDCProc,'EfficientSubpixelRegistrationProcess')
        meanYShift = round(T(1,1));
        meanXShift = round(T(1,2));
        firstMask = circshift(tempMask,[meanYShift meanXShift]);
        % Now I blacked out erroneously circularaly shifted bead image
        % portion - SH 20171008
        if meanYShift>=0 %shifted downward
            firstMask(1:meanYShift,:)=0;
        else %shifted upward
            firstMask(end+meanYShift:end,:)=0;
        end
        if meanXShift>=0 %shifted right hand side
            firstMask(:,1:meanXShift)=0;
        else %shifted left
            firstMask(:,end+meanXShift:end)=0;
        end
    else
        %y_shift = find(any(firstMask,2),1);
        %x_shift = find(any(firstMask,1),1);
        %tempMask2(y_shift:y_shift+size(tempMask,1)-1,x_shift:x_shift+size(tempMask,2)-1) = tempMask;
        Tr = affine2d([1 0 0; 0 1 0; fliplr(T(1, :)) 1]);
        tempMask2 = imwarp(tempMask,Tr,'OutputView',ref_obj);
        firstMask = tempMask2 & firstMask;
    end
end
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
debuggingMode = false; %true; % default should be false. If you want use it, change currentBeads at line 349
for j= firstFrame:nFrames
    % Read image and perform correlation
    if ~isempty(SDCProc)
        currImage = double(SDCProc.loadChannelOutput(p.ChannelIndex(1),j));
    else
        currImage = double(movieData.channels_(p.ChannelIndex(1)).loadImage(j));
    end
    if ~p.useGrid
    % if strcmp(movieData.getChannel(p.ChannelIndex).imageType_,'Widefield')
        % Detect beads in reference frame
        if j==firstFrame && firstFrame==1
            disp('Determining PSF sigma from reference frame...')
            if ~debuggingMode
                % psfSigma = getGaussianSmallestPSFsigmaFromData(refFrame,'Display',false);
%                 psfSigma=NaN;
%                 if isnan(psfSigma) || psfSigma>movieData.channels_(p.ChannelIndex(1)).psfSigma_*3 
%                     if strcmp(movieData.getChannel(p.ChannelIndex(1)).imageType_,'Widefield') || movieData.pixelSize_>130
%                         psfSigma = movieData.channels_(p.ChannelIndex(1)).psfSigma_*2; %*2 scale up for widefield
%                     elseif strcmp(movieData.getChannel(p.ChannelIndex(1)).imageType_,'Confocal')
%                         psfSigma = movieData.channels_(p.ChannelIndex(1)).psfSigma_*0.79; %*4/7 scale down for  Confocal finer detection SH012913
%                     elseif strcmp(movieData.getChannel(p.ChannelIndex(1)).imageType_,'TIRF')
%                         psfSigma = movieData.channels_(p.ChannelIndex(1)).psfSigma_*4/7; %*3/7 scale down for TIRF finer detection SH012913
%                     else
%                         disp('image type should be chosen among Widefield, confocal and TIRF! Assuming TIRF for now...');
%                         psfSigma = movieData.channels_(p.ChannelIndex(1)).psfSigma_*4/7; %*3/7 scale down for TIRF finer detection SH012913
%                     end
%                 end
                psfSigma = movieData.channels_(p.ChannelIndex(1)).psfSigma_;
                disp(['Determined sigma: ' num2str(psfSigma)])

                disp('Detecting beads in the reference frame...')
                try
                    pstruct = pointSourceDetection(refFrame, psfSigma, 'alpha', p.alpha,'Mask',firstMask); %,...
                        %'FitMixtures',true, 'MaxMixtures', 3);
                catch
                    % in this case, firstMask is wrong due to cropped
                    % image. Setting it again...
                    if ~isempty(SDCProc)
                        firstImage = double(SDCProc.loadChannelOutput(p.ChannelIndex(1),1));
                    else
                        firstImage = double(movieData.channels_(p.ChannelIndex(1)).loadImage(1));
                    end
                    firstMask = refFrame>0 & firstImage>0;
                    pstruct = pointSourceDetection(refFrame, psfSigma, 'alpha', p.alpha,'Mask',firstMask,...
                        'FitMixtures',true, 'MaxMixtures', 3);
                    
                end
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
                    avg_beads_distance = quantile(distance,0.5);%mean(distance);%size(refFrame,1)*size(refFrame,2)/length(beads);
                    notSaturated = true;
                    xmin = min(pstruct.x);
                    xmax = max(pstruct.x);
                    ymin = min(pstruct.y);
                    ymax = max(pstruct.y);
                %     avgAmp = mean(pstruct.A);
                %     avgBgd = mean(pstruct.c);
                %     thresInten = avgBgd+0.02*avgAmp;
    %                 thresInten = quantile(pstruct.c,0.25);
                    thresInten = quantile(pstruct.c+pstruct.A,0.5); % try to pick up bright-enough spots
                    maxNumNotDetected = 20; % the number of maximum trial without detecting possible point
                    numNotDetected = 0;
                    numPrevBeads = size(beads,1);
                    % To avoid camera noise, Gaussian-filtered image will be
                    % used - SH 20171010
                    refFrameFiltered = filterGauss2D(refFrame,psfSigma*0.9);
                    tic
                    while notSaturated
                        x_new = xmin + (xmax-xmin)*rand(10000,1);
                        y_new = ymin + (ymax-ymin)*rand(10000,1);
                        [~,distToPoints] = KDTreeClosestPoint(beads,[x_new,y_new]);
                        inten_new = arrayfun(@(x,y) refFrameFiltered(round(y),round(x)),x_new,y_new);
                        idWorthAdding = distToPoints>avg_beads_distance/3 & inten_new>thresInten;
                        if sum(idWorthAdding)>1
                            beads = [beads; [x_new(idWorthAdding), y_new(idWorthAdding)]];
                            numNotDetected = 0;
                        else
                            numNotDetected=numNotDetected+1;
                        end
                        if numNotDetected>maxNumNotDetected
                            notSaturated = false; % this means now we have all points to start tracking from the image
                        end
                    end
                    toc
                    disp([num2str(size(beads,1)-numPrevBeads) ' points were additionally detected for fine tracking. Total detected beads: ' num2str(length(beads))])
                end
                % Exclude all beads which are less  than half the correlation length 
                % away from the padded border. By default, no centered template should 
                % include any NaN's for correlation
                % Create beads mask with zero intensity points as false
                beadsMask = firstMask; %true(size(refFrame));
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
                beadsMask2=beadsMask & bwmorph(true(size(refFrame)),'erode',erosionDist);
                %         beadsMask=imerode(beadsMask,strel('square',erosionDist));
                indx=beadsMask2(sub2ind(size(beadsMask),ceil(beads(:,2)), ceil(beads(:,1))));
                localbeads = beads(indx,:);
                currentBeads = localbeads; %This will keep updated
                cumulativeV_forV=zeros(size(localbeads));
                cumulativeV_forBeads=zeros(size(localbeads));
            else % Debugging mode
                localbeads = [545 880];
                currentBeads = localbeads; %This will keep updated
            end
        elseif j==firstFrame && firstFrame>1
            localbeads = displField(j-1).pos;
            currentBeads = localbeads; %This will keep updated
        end
        %     % Select only beads which are min correlation length away from the border of the
        %     % reference frame 
        %     beadsMask = true(size(refFrame));
        %     erosionDist=p.minCorLength;
        %     % erosionDist=p.minCorLength+1;
        %     beadsMask(erosionDist:end-erosionDist,erosionDist:end-erosionDist)=false;
        %     indx=beadsMask(sub2ind(size(beadsMask),ceil(beads(:,2)),ceil(beads(:,1))));
        %     beads(indx,:)=[];
        if ~p.trackSuccessively
            scoreCalculation='xcorr';
            % Track beads displacement in the xy coordinate system
            v = trackStackFlow(cat(3,refFrame,currImage),currentBeads,...
                p.minCorLength,p.minCorLength,'maxSpd',p.maxFlowSpeed,...
                'mode',p.mode,'scoreCalculation',scoreCalculation,'sigCrit',p.sigCrit,'noGradualExpansionOfSearchArea',p.noGradualExpansionOfSearchArea>0);%,'usePIVSuite', p.usePIVSuite);
        else
%             scoreCalculation='difference';
            scoreCalculation='xcorr';
            % Track beads displacement in the xy coordinate system
            if j==firstFrame
                prevImage=refFrame;
            end
            v = trackStackFlow(cat(3,prevImage,currImage),currentBeads,...
                p.minCorLength,p.minCorLength,'maxSpd',p.maxFlowSpeed,...
                'mode',p.mode,'scoreCalculation',scoreCalculation,'sigCrit',p.sigCrit);
            prevImage=currImage;
        end

        % Extract finite displacement and prepare displField structure in the xy
        % coordinate system
%         if j==firstFrame
        validV = ~isinf(v(:,1)) & ~isnan(v(:,1));
%         end
        
        if ~p.trackSuccessively
%             displField(j).pos=localbeads(validV,:);
%             displField(j).vec=[v(validV,1)+residualT(j,2) v(validV,2)+residualT(j,1)]; % residual should be added with oppiste order! -SH 072514
            displField(j).pos=localbeads; % validV is removed to include NaN location - SH 030417
            displField(j).vec=[v(:,1)+residualT(j,2) v(:,2)+residualT(j,1)]; % residual should be added with oppiste order! -SH 072514
        else
            v2 = v;
            cumulativeV_forV = cumulativeV_forV+v2;
            v2(~validV,1)=0; v2(~validV,2)=0;
            cumulativeV_forBeads = cumulativeV_forBeads+v2;
            currentBeads = localbeads + cumulativeV_forBeads;
            displField(j).pos=localbeads(validV,:);
            displField(j).vec=[cumulativeV_forV(validV,1)+residualT(j,2) cumulativeV_forV(validV,2)+residualT(j,1)]; % residual should be added with oppiste order! -SH 072514
        end
    else
        % Just in case there is a path problem
        addpath('~/matlab/extern/PIVsuite/PIVsuite v.0.81')
        
        pivPar = [];      % variable for settings
        pivData = [];     % variable for storing results

        [pivPar, pivData] = pivParams(pivData,pivPar,'defaults');     
        % Set the size of interrogation areas via fields |iaSizeX| and |iaSizeY| of |pivPar| variable:
%         pivPar.iaSizeX = [64 32 16 2^(nextpow2(p.minCorLength)-1)];     % size of interrogation area in X 
        nextPow2max=nextpow2(p.maxFlowSpeed);
        nextPow2corLeng=nextpow2(p.minCorLength);
        
        sizeArray=2.^([nextPow2max nextPow2max nextPow2max nextPow2max]);
        stepArray=2.^([nextPow2corLeng nextPow2corLeng:-1:nextPow2corLeng-2]);
        pivPar.anNpasses = length(stepArray);
        pivPar.iaSizeX = sizeArray;     % size of interrogation area in X 
        pivPar.iaStepX = stepArray;     % grid spacing of velocity vectors in X
        pivPar.iaSizeY = sizeArray;     % size of interrogation area in X 
        pivPar.iaStepY = stepArray;    % grid spacing of velocity vectors in X
%         pivPar.iaSizeX = [64 32 16 8];     % size of interrogation area in X 
%         pivPar.iaStepX = [32 16  8 4];     % grid spacing of velocity vectors in X
%         pivPar.iaSizeY = [64 32 16 8];     % size of interrogation area in X 
%         pivPar.iaStepY = [32 16  8 4];     % grid spacing of velocity vectors in X
% %         pivPar.iaStepX = [32 16  8 8];     % grid spacing of velocity vectors in X
% %         pivPar.iaSizeY = [64 32 16 16];     % size of interrogation area in Y 
% %         pivPar.iaStepY = [32 16  8 8];     % grid spacing of velocity vectors in Y
        pivPar.ccMaxDisplacement = 0.7;%p.maxFlowSpeed/sizeArray(end);   % This filter is relatively narrow and will 
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
%         validV = ~isnan(pivData.V);
        validV = true(size(pivData.V));
        
        if ~p.trackSuccessively
            displField(j).pos=[pivData.X(validV), pivData.Y(validV)];
            displField(j).vec=[pivData.U(validV)+residualT(j,2), pivData.V(validV)+residualT(j,1)]; % residual should be added with oppiste order! -SH 072514
        else
            v2 = [pivData.U, pivData.V];
            if j== firstFrame
                cumulativeV_forV = zeros(size(pivData.U,1),2);
                cumulativeV_forBeads = zeros(size(pivData.U,1),2);
            else
                cumulativeV_forV = cumulativeV_forV+v2;
                v2(~validV,1)=0; v2(~validV,2)=0;
                cumulativeV_forBeads = cumulativeV_forBeads+v2;
            end
            currentBeads = [pivData.X(validV), pivData.Y(validV)] + cumulativeV_forBeads;
            displField(j).pos=[pivData.X(validV), pivData.Y(validV)];
            displField(j).vec=[cumulativeV_forV(validV,1)+residualT(j,2) cumulativeV_forV(validV,2)+residualT(j,1)]; % residual should be added with oppiste order! -SH 072514
        end

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
    if ishandle(wtBar) %mod(j,5)==1 && 
        tj=toc;
        waitbar(j/nFrames,wtBar,sprintf([logMsg ...
            timeMsg(tj*(nFrames-firstFrame+1-j)/j)]));
    end
    disp(['Done for frame # ' num2str(j)])
    
    % Save each iteration (for recovery of unfinished processes)
    save(outputFile{1},'displField');
end
%% Displacement map creation - this is shifted version
% Now decided to discard this part to save the hard disk space 20180421 SH
% [dMapIn, dmax, dmin, cropInfo,dMapXin,dMapYin,reg_grid] = generateHeatmapShifted(displField,displField,0);
% disp(['Estimated displacement maximum = ' num2str(dmax) ' pixel.'])
% Insert traction map in forceField.pos 
% disp('Generating displacement maps ...')
% outputDir = fullfile(p.OutputDirectory,'displMaps');
% mkClrDir(outputDir);
% fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];
% numStr = @(frame) num2str(frame,fString);
% outFileDMap=@(frame) [outputDir filesep 'displMap' numStr(frame) '.mat'];

% dMap = cell(1,nFrames);
% dMapX = cell(1,nFrames);
% dMapY = cell(1,nFrames);
[reg_grid,~,~,~]=createRegGridFromDisplField(displField,1.0,0);
displFieldShifted(nFrames)=struct('pos','','vec','');
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
save(outputFile{1},'displField','displFieldShifted','-v7.3');% I removed v7.3 option to save the space,'-v7.3');
% Saving the dMap which stores information
% dMap.eachDMapName = 'cur_dMap';
% dMap.outputDir = fullfile(p.OutputDirectory,'dislplMaps');
% dMap.outFileDMap = @(frame) [outputDir filesep 'tractionMap' numStr(frame) '.mat'];
dMap.displFieldPath = outputFile{1};
dMap.firstMaskSize = size(firstMask);

dMapX=dMap; %dMapX.eachTMapName = 'cur_dMapX';
dMapY=dMap; %dMapY.eachTMapName = 'cur_dMapY';
save(outputFile{2},'dMap','dMapX','dMapY'); % Updated, SH 20180225
% save(outputFile{2},'dMap','dMapX','dMapY','-v7.3'); % need to be updated for faster loading. SH 20141106
displMag=cell2mat(arrayfun(@(x) sqrt(x.vec(:,1).^2+x.vec(:,2).^2), displField,'unif',false)');
dmin = quantile(displMag(:),0.01); dmax = quantile(displMag(:),0.95);
displFieldProc.setTractionMapLimits([dmin dmax])

% Close waitbar
if ishandle(wtBar), close(wtBar); end

disp('Finished calculating displacement field!')