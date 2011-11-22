function trackMovieFlow(movieData,varargin)
% trackMovieFlow track the flow of a movie using image correlation
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

% Sebastien Besson 5/2011 (last modified Nov 2011)

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
else
    wtBar=-1;
end

%Find the speckle detection process, and the segmentation process
iMaskProc =movieData.getProcessIndex('MaskRefinementProcess',1,1);
iSpecProc =movieData.getProcessIndex('SpeckleDetectionProcess',1,1);

nChan = length(p.ChannelIndex);
if isempty(iSpecProc) || isempty(iMaskProc)
    error(['Speckle detection and segmentation have not yet been performed '...
        'on this movie! Please run first!!']);
end
segProc = movieData.processes_{iMaskProc};
specDetProc = movieData.processes_{iSpecProc};

%Check which channels have speckles and masks
hasMasks = segProc.checkChannelOutput(p.ChannelIndex);
hasSpec = specDetProc.checkChannelOutput(p.ChannelIndex);
if ~all(hasMasks & hasSpec)
    error(['Each channel must have speckles and masks! ' ...
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
firstFrameIndx = p.firstImage:p.timeStepSize:p.lastImage -p.timeWindow+1;
nFrames = numel(firstFrameIndx);
nTot = nChan*nFrames;

logMsg = @(chan) ['Please wait, tracking flow for channel ' num2str(chan)];
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t/60)) 'min'];
tic;

for i = 1:numel(p.ChannelIndex)
    iChan = p.ChannelIndex(i);
    % Log display
    disp(logMsg(iChan))
    disp(movieData.getChannelPaths{iChan});
    disp('Results will be saved under:')
    disp(outFilePaths{1,iChan});
    
    % Load raw images and masks from segmentation output
    if ishandle(wtBar), waitbar(0,wtBar,'Loading images and masks...'); end
    disp('Loading images and masks...');
    stack = zeros([movieData.imSize_ movieData.nFrames_]);
    bgMask = zeros([movieData.imSize_ movieData.nFrames_]);
    maskNames = segProc.getOutMaskFileNames(iChan);
    inMask=@(frame) [flowTrackProc.inFilePaths_{2,iChan}...
        filesep maskNames{1}{frame}];
    for j = p.firstImage:p.lastImage
        stack(:,:,j) = movieData.channels_(iChan).loadImage(j);
        bgMask(:,:,j) = imerode(logical(imread(inMask(j))),...
            strel('disk',p.edgeErodeWidth));
    end
    
    % Load speckles from speckle detection output
    if ishandle(wtBar), waitbar(0,wtBar,'Loading speckles...'); end
    disp('Loading speckles...');
    speckles = cell(1,nFrames);
    for j = firstFrameIndx
        cands = specDetProc.loadChannelOutput(iChan,j);
        speckles{j} = vertcat(cands([cands.status]==1).Lmax);
    end
    
    % Call the main correlation routine
    disp('Starting tracking flow...');
    if ishandle(wtBar), waitbar(0,wtBar,logMsg(iChan)); end
    for j=firstFrameIndx
        % Call the flow tracking routine
        [vx,vy,corLen] = trackStackFlow(stack(:,:,j:j+p.timeWindow-1),...
            speckles{j}(:,2),speckles{j}(:,1),p.minCorLength,p.maxCorLength,...
            'maxSpd',p.maxFlowSpeed,'bgMask',bgMask(:,:,j:j+p.timeWindow-1), ...
            'numStBgForAvg',p.numStBgForAvg,'minFeatureSize',p.minFeatureSize);  %#ok<NASGU>
        
        % Concatenate flow as a [pos1 pos2] matrix
        flow = [speckles{j} vy vx];
        corLen = [speckles{j} corLen];
        % Set infinite flow to nan
        flow(isinf(vx),3:4)=NaN;
                
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
        if mod(j,5)==1 && ishandle(wtBar), 
            tj=toc;
            nj = (i-1)*nFrames+ j;
            waitbar(nj/nTot,wtBar,sprintf([logMsg(iChan) timeMsg(tj*nTot/nj-tj)]));
        end
    end
end

% Close waitbar
if ishandle(wtBar), close(wtBar); end
disp('Finished tracking flow!')