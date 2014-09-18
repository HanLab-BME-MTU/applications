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
outputFile{1} = [p.OutputDirectory filesep 'displField.mat'];

% Add a recovery mechanism if process has been stopped in the middle of the
% computation to re-use previous results
firstFrame =1; % Set the strating fram eto 1 by default
if exist(outputFile{1},'file');
    % Check analyzed frames
    s=load(outputFile{1},'displField');
    frameDisplField=~arrayfun(@(x)isempty(x.pos),s.displField);
    
    if ~all(frameDisplField) && ~all(~frameDisplField)
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
% Detect beads in reference frame 
disp('Detecting beads in the reference frame...')
if strcmp(movieData.getChannel(p.ChannelIndex).imageType_,'Widefield')
    sigmaPSF = movieData.channels_(1).psfSigma_*2; %*2 scale up for widefield
elseif strcmp(movieData.getChannel(p.ChannelIndex).imageType_,'Confocal')
    sigmaPSF = movieData.channels_(1).psfSigma_*0.79; %*4/7 scale down for  Confocal finer detection SH012913
elseif strcmp(movieData.getChannel(p.ChannelIndex).imageType_,'TIRF')
    sigmaPSF = movieData.channels_(1).psfSigma_*3/7; %*3/7 scale down for TIRF finer detection SH012913
else
    error('image type should be chosen among Widefield, confocal and TIRF!');
end
pstruct = pointSourceDetection(refFrame, sigmaPSF, 'alpha', p.alpha,'Mask',firstMask);
assert(~isempty(pstruct), 'Could not detect any bead in the reference frame');
beads = [ceil(pstruct.x') ceil(pstruct.y')];

% Subsample detected beads ensuring beads are separated by at least half of
% the correlation length - commented out to get more beads
if ~p.highRes
    disp('Subsampling detected beads (normal resolution)...')
    max_beads_distance = floor(p.minCorLength/2);
else
    % To get high-resolution information, subsample detected beads ensuring 
    % beads are separated by 0.1 um the correlation length 
    disp('Subsampling detected beads (high resolution)...')
    max_beads_distance = floor(100/movieData.pixelSize_);
end

idx = KDTreeBallQuery(beads, beads, max_beads_distance);
valid = true(numel(idx),1);
for i = 1 : numel(idx)
    if ~valid(i), continue; end
    neighbors = idx{i}(idx{i} ~= i);
    valid(neighbors) = false;
end
beads = beads(valid, :);

% Select only beads which are min correlation length away from the border of the
% reference frame 
beadsMask = true(size(refFrame));
erosionDist=p.minCorLength;
% erosionDist=p.minCorLength+1;
beadsMask(erosionDist:end-erosionDist,erosionDist:end-erosionDist)=false;
indx=beadsMask(sub2ind(size(beadsMask),ceil(beads(:,2)),ceil(beads(:,1))));
beads(indx,:)=[];

if p.useGrid && p.highRes
    tempDisplField.pos = beads;
    [~,xvec,yvec,~]=createRegGridFromDisplField(tempDisplField,2,0);
    beads = [xvec yvec];
end

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
    % Exclude all beads which are less  than half the correlation length 
    % away from the padded border. By default, no centered template should 
    % include any NaN's for correlation
    % Create beads mask with zero intensity points as false
    beadsMask = true(size(refFrame));
    beadsMask(currImage==0)=false;
    % Remove false regions non-adjacent to the image border
    beadsMask = beadsMask | imclearborder(~beadsMask);
    % Erode the mask with half the correlation length and filter beads
    erosionDist=round((p.minCorLength+1)/2);
    beadsMask=imerode(beadsMask,strel('square',erosionDist));
    indx=beadsMask(sub2ind(size(beadsMask),ceil(beads(:,2)), ceil(beads(:,1))));
    localbeads = beads(indx,:);

    % Track beads displacement in the xy coordinate system
    v = trackStackFlow(cat(3,refFrame,currImage),localbeads,...
        p.minCorLength,p.minCorLength,'maxSpd',p.maxFlowSpeed,...
        'mode',p.mode);
    
    % Extract finite displacement and prepare displField structure in the xy
    % coordinate system
    validV = ~isinf(v(:,1));
    displField(j).pos=localbeads(validV,:);
    displField(j).vec=[v(validV,1)+residualT(j,2) v(validV,2)+residualT(j,1)]; % residual should be added with oppiste order! -SH 072514
    
    % Update the waitbar
    if mod(j,5)==1 && ishandle(wtBar)
        tj=toc;
        waitbar(j/nFrames,wtBar,sprintf([logMsg ...
            timeMsg(tj*(nFrames-firstFrame++1-j)/j)]));
    end
    
    % Save each iteration (for recovery of unfinished processes)
    save(outputFile{1},'displField');
end

% Close waitbar
if ishandle(wtBar), close(wtBar); end

disp('Finished calculating displacement field!')