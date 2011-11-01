function trackMovieFlow(movieData,varargin)
% trackMovieFlow track the flow of a movie using image correlation
%
%
% SYNOPSIS trackMovieFlow(movieData,paramsIn)
%
% INPUT
%   movieData - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%

% Sebastien Besson 5/2011

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous speckle detection processes
iProc = movieData.getProcessIndex('FlowTrackingProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(FlowTrackingProcess(movieData,...
        movieData.outputDirectory_));
end
flowTrackProc = movieData.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(flowTrackProc,paramsIn);

%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...');
end

%Find the speckle detection process, and the segmentation process
iSegProc =movieData.getProcessIndex('MaskRefinementProcess',1,1);
iSpecProc =movieData.getProcessIndex('SpeckleDetectionProcess',1,1);

nChan = length(p.ChannelIndex);
if isempty(iSpecProc) || isempty(iSegProc)
    error(['Speckle detection and segmentation have not yet been performed '...
        'on this movie! Please run first!!']);
end
segProc = movieData.processes_{iSegProc};
specDetProc = movieData.processes_{iSegProc};

%Check which channels have speckles and masks
hasMasks = segProc.checkChannelOutput(p.ChannelIndex);
hasSpeck = specDetProc.checkChannelOutput(p.ChannelIndex);
if ~all(hasMasks && hasSpeck)
    error(['Each channel must have speckles! ' ...
        'Please apply speckle detection to all needed channels before '...
        'running flow tracking!'])
end


% Set up the input directories
inFilePaths = cell(1,1:numel(movieData.channels_));
for i = p.ChannelIndex
    inFilePaths{1,i} = specDetProc.outFilePaths_{1,i};
    inFilePaths{2,i} = segProc.outFilePaths_{1,i};
end
flowTrackProc.setInFilePaths(inFilePaths);

% Set up the output directories
outFilePaths = cell(2,1:numel(movieData.channels_));
for i = p.ChannelIndex;
    %Create string for current directory
    outFilePaths{1,i} = [p.OutputDirectory filesep 'flow_for_channel_' num2str(i)];
    mkClrDir(outFilePaths{1,i});
end
flowTrackProc.setOutFilePaths(outFilePaths);

%% --------------- Flow tracking ---------------%%%

disp('Starting tracking flow...')
% Reading various constants
imDirs  = movieData.getChannelPaths;
imageFileNames = movieData.getImageFileNames;
nFrames = movieData.nFrames_;

% Anonymous functions for reading input/output
inImage=@(chan,frame) [imDirs{chan} filesep imageFileNames{chan}{frame}];

logMsg = @(chan) ['Please wait, tracking flow for channel ' num2str(chan)];
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t)) 's'];
tic;
nTot = nChan*nFrames;

for i = 1:numel(p.ChannelIndex)
    iChan = p.ChannelIndex(i);
    % Log display
    disp(logMsg(iChan))
    disp(imDirs{iChan});
    disp('Results will be saved under:')
    disp(outFilePaths{1,iChan});
    
    % Load raw images and masks from segmentation output
    if feature('ShowFigureWindows')
        waitbar(0,wtBar,'Loading images and masks...');
    end
    disp('Loading images and masks...');
    stack = zeros([movieData.imSize_ movieData.nFrames_]);
    bgMask = zeros([movieData.imSize_ movieData.nFrames_]);
    maskNames = segProc.getOutMaskFileNames(iChan);
    inMask=@(frame) [flowTrackProc.inFilePaths_{2,iChan}...
        filesep maskNames{1}{frame}];
    for j = p.firstImage:p.lastImage
        stack(:,:,j) = imread(inImage(iChan,j));
        bgMask(:,:,j) = imerode(logical(imread(inMask(j))),...
            strel('disk',p.edgeErodeWidth));
    end
    
    % Load speckles from speckle detection output
    if feature('ShowFigureWindows')
        waitbar(0,wtBar,'Loading speckles...');
    end
    disp('Loading speckles...');
    firstImagesIndx = p.firstImage:p.timeStepSize:p.lastImage -p.timeWindow+1;
    speckles = cell(1,length(firstImagesIndx));
    for j = firstImagesIndx
        cands = movieData.processes_{iSpecProc}.loadChannelOutput(iChan,j);
        M = vertcat(cands([cands.status]==1).Lmax);
        speckles{j} = M(:,2:-1:1);
    end
    
    % Call the main correlation routine
    disp('Starting tracking flow...');
    if feature('ShowFigureWindows'), waitbar(0,wtBar,logMsg(iChan)); end
    for j=firstImagesIndx
        [vx,vy,corLen] = trackStackFlow(stack(:,:,j:j+p.timeWindow-1),...
            speckles{j}(:,1),speckles{j}(:,2),...
            p.minCorLength,p.maxCorLength,...
            'maxSpd',p.maxFlowSpeed,'bgMask',bgMask(:,:,j:j+p.timeWindow-1), ...
            'numStBgForAvg',p.numStBgForAvg,'minFeatureSize',p.minFeatureSize);
        velocity = [vx vy];
        
        % Concatenate flow as a [pos1 pos2] matrix
        flow = [speckles{j}(:,2:-1:1) velocity(:,2:-1:1,1)];
        
        % Set infinite flow to nan
        finiteFlow  = ~isinf(velocity(:,1));
        flow(~finiteFlow,3:4)=NaN;
        
        % Filter vector field outliers
        if ~isempty(p.outlierThreshold)
            outlierIndex = detectVectorFieldOutliers(flow,p.outlierThreshold);
            flow(outlierIndex,3:4)=NaN;
        end
        
        % Save flow result under [pos pos+vel] format
        flow(:,3:4)=flow(:,1:2)+flow(:,3:4); %#ok<NASGU>
        flowFileName = [outFilePaths{1,iChan} filesep 'flow' num2str(j) '.mat'];
        save(flowFileName,'flow','corLen');
        
        % Update waitbar
        if mod(j,5)==1 && feature('ShowFigureWindows')
            tj=toc;
            nj = (i-1)*nFrames+ j;
            waitbar(nj/nTot,wtBar,sprintf([logMsg(iChan) timeMsg(tj*nTot/nj-tj)]));
        end

    end
end
% Close waitbar
if feature('ShowFigureWindows'), close(wtBar);end
disp('Finished tracking flow!')