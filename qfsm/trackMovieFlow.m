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
wtBar=-1;
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
maskProc = movieData.processes_{iMaskProc};
specDetProc = movieData.processes_{iSpecProc};

%Check which channels have speckles and masks
hasMasks = maskProc.checkChannelOutput(p.ChannelIndex);
hasSpec = specDetProc.checkChannelOutput(p.ChannelIndex);
if ~all(hasMasks & hasSpec)
    error(['Each channel must have speckles and masks! ' ...
        'Please apply speckle detection to all needed channels before '...
        'running flow tracking!'])
end


% Set up the input directories
inFilePaths = cell(3,numel(movieData.channels_));
imDirs = movieData.getChannelPaths;
for i = p.ChannelIndex
    inFilePaths{1,i} = imDirs{i};
    inFilePaths{2,i} = specDetProc.outFilePaths_{1,i};
    inFilePaths{3,i} = maskProc.outFilePaths_{1,i};
end
flowTrackProc.setInFilePaths(inFilePaths);

% Set up the output directories
outFilePaths = cell(2,numel(movieData.channels_));
channelName = @(x)movieData.getChannelPaths{x}(max(regexp(movieData.getChannelPaths{x},filesep))+1:end);
for i = p.ChannelIndex;
    %Create string for current directory
    outFilePaths{1,i} = fullfile(p.OutputDirectory,channelName(i));
    mkClrDir(outFilePaths{1,i});
end
flowTrackProc.setOutFilePaths(outFilePaths);

%% --------------- Flow tracking ---------------%%%

disp('Starting tracking flow...')
% Reading various constants
if p.timeStepSize==0
    firstFrames = p.firstImage;
else
    firstFrames = p.firstImage:p.timeStepSize:p.lastImage-p.timeWindow+1;
end 
imSize= movieData.imSize_;
nFrames = movieData.nFrames_;
nStacks=numel(firstFrames);
nTot = nChan*nStacks;

% Get number of stationary images for background
if p.numStBgForAvg==-1, 
    p.numStBgForAvg=nFrames; 
elseif p.numStBgForAvg>0,
    p.numStBgForAvg = min(nFrames,p.numStBgForAvg);
end

logMsg = @(chan) ['Please wait, tracking flow for channel ' num2str(chan)];
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t/60)) 'min'];
tic;

channelLog=cell(numel(p.ChannelIndex),1);
for i = 1:numel(p.ChannelIndex)
    iChan = p.ChannelIndex(i);
     % Log display
    channelLog{i} = sprintf('Channel %g: %s\n',iChan,imDirs{iChan});
    disp(logMsg(iChan))
    disp('Using images from directory:')
    disp(inFilePaths{1,iChan});
    disp('and speckles from diretory:')
    disp(inFilePaths{2,iChan});
    disp('Results will be saved under:')
    disp(outFilePaths{1,iChan});
    
    % Initialize input
    stack = zeros([imSize nFrames]);
    bgMask = true([imSize nFrames]);
    bgAvgImg=zeros([imSize nFrames]);
    
    % Load raw images (incl. images for background averaging)
    if ishandle(wtBar), waitbar(0,wtBar,'Loading images ...'); end
    disp('Loading images...');
    minFrame=min(p.firstImage,nFrames-p.numStBgForAvg+1);
    maxFrame=max(p.lastImage,min(p.lastImage-p.timeWindow+1+p.numStBgForAvg,nFrames));
    for j = minFrame:maxFrame
        stack(:,:,j) = movieData.channels_(iChan).loadImage(j);
    end
    
    % Construct background stack by averaging stationary images
    if p.numStBgForAvg>0
        nStBgFrames=p.timeWindow-min(p.timeWindow,p.numStBgForAvg)+1;
        for j = firstFrames
            for k=j:j+nStBgFrames-1
                startStBgFrame = min(k,nFrames-p.numStBgForAvg+1);
                bgAvgImg(:,:,k)=mean(stack(:,:,startStBgFrame:startStBgFrame+p.numStBgForAvg-1),3);
            end
        end
    else
        nStBgFrames=1;
    end

    
    % Load masks
    if ishandle(wtBar), waitbar(0,wtBar,'Loading masks ...'); end
    disp('Loading masks...');
    se=strel('disk',p.edgeErodeWidth);
    for j = p.firstImage:p.lastImage
        bgMask(:,:,j) = imerode(maskProc.loadChannelOutput(iChan,j),se);
    end
    
    if ishandle(wtBar), waitbar(0,wtBar,'Loading speckles...'); end
    disp('Loading speckles...');
    speckles = cell(1,nFrames);
    for firstFrame = firstFrames
        % Load speckles from speckle detection output
        cands = specDetProc.loadChannelOutput(iChan,firstFrame);
        % Retrieve position of significant candidates (in image coordinate system)
        speckles{firstFrame} = vertcat(cands([cands.status]==1).Lmax);
    end
    
    % Call the main correlation routine
    disp('Starting tracking flow...');
    if ishandle(wtBar), waitbar(0,wtBar,logMsg(iChan)); end
    allCorLen=cell(numel(firstFrames));
    
    for j=1:numel(firstFrames)
        firstFrame=firstFrames(j);
        frameRange =firstFrame:firstFrame+p.timeWindow-1;
        
        if p.usePIV
            % Just in case there is a path problem
            addpath('~/matlab/extern/PIVsuite/PIVsuite v.0.81')

            pivPar = [];      % variable for settings
            pivData = [];     % variable for storing results

            [pivPar, pivData] = pivParams(pivData,pivPar,'defaults');     
            % Set the size of interrogation areas via fields |iaSizeX| and |iaSizeY| of |pivPar| variable:
    %         pivPar.iaSizeX = [64 32 16 2^(nextpow2(p.minCorLength)-1)];     % size of interrogation area in X 
            nextPow2max=nextpow2(p.maxCorLength);
            nextPow2min=nextpow2(p.minCorLength);

            sizeArray=2.^([nextPow2max nextPow2max round((nextPow2max+nextPow2min)/2) nextPow2min]);
            stepArray=sizeArray;
            pivPar.anNpasses = length(stepArray);
            pivPar.iaSizeX = sizeArray;     % size of interrogation area in X 
            pivPar.iaStepX = stepArray;     % grid spacing of velocity vectors in X
            pivPar.iaSizeY = sizeArray;     % size of interrogation area in X 
            pivPar.iaStepY = stepArray;    % grid spacing of velocity vectors in X
            pivPar.ccMaxDisplacement = 0.7;%p.maxFlowSpeed/sizeArray(end);   % This filter is relatively narrow and will 
            %pivPar.ccWindow = 'Gauss2';   %will use Welch, which is default.
            pivPar.smMethod = 'none';
            pivPar.iaMethod = 'defspline';
            pivPar.iaImageInterpolationMethod = 'spline';
            pivPar.imMask1=bgMask(:,:,frameRange(1));
            pivPar.imMask2=bgMask(:,:,frameRange(end));

            [pivData] = pivAnalyzeImagePair(stack(:,:,frameRange(1)),...
                stack(:,:,frameRange(end)),pivData,pivPar);
            validV = true(size(pivData.V));

            flow=[pivData.Y(validV) pivData.X(validV) ...
                pivData.V(validV) pivData.U(validV)];       
            
            corLen = pivData.iaSizeX;
            allCorLen{j}=corLen;
        else
            % Call the main flow tracking routine using speckles positions in xy
            % coordinate system
            [v,corLen] = trackStackFlow(stack(:,:,frameRange),...
                speckles{firstFrame}(:,2:-1:1),p.minCorLength,p.maxCorLength,...
                'maxSpd',p.maxFlowSpeed,'bgMask',bgMask(:,:,frameRange), ...
                'bgAvgImg',bgAvgImg(:,:,firstFrame:firstFrame+nStBgFrames-1),...
                'minFeatureSize',p.minFeatureSize);
            allCorLen{j}=corLen;
            % Concatenate flow as a [pos1 pos2] matrix into image coordinate system
            flow = [speckles{firstFrame} v(:,2:-1:1)];

            % Set infinite flow to nan
            flow(isinf(v(:,1)),3:4)=NaN;
        end        
                
        % Filter vector field outliers
        if ~isempty(p.outlierThreshold)
            outlierIndex = detectVectorFieldOutliers(flow,p.outlierThreshold);
            flow(outlierIndex,3:4)=NaN;
        end
        
        % Save flow result under [pos pos+vel] format
        flow(:,3:4)=flow(:,1:2)+flow(:,3:4); %#ok<NASGU>
        flowFileName = [outFilePaths{1,iChan} filesep 'flow' num2str(firstFrame) '.mat'];
        save(flowFileName,'flow','corLen');
        
        % Update waitbar
        if mod(j,5)==1 && ishandle(wtBar), 
            tj=toc;
            nj = (i-1)*nStacks+ j;
            waitbar(nj/nTot,wtBar,sprintf([logMsg(iChan) timeMsg(tj*nTot/nj-tj)]));
        end
    end
    
    % Get count of finit correlation lengths
    allCorLen=vertcat(allCorLen{:});
    finiteCorLen =allCorLen(~isinf(allCorLen));
    [corLenValues,~,index]=unique(finiteCorLen);
    corLenCount = hist(index,unique(index));
    
    % Create channel log fot output
    channelLog{i} = [channelLog{i} ...
        sprintf([...
        'Integration window\t\t\t\t\t\t: %g\n'...
        'Total number of speckles in the movie\t\t\t\t: %g\n'...
        'Total number of tracked points overall\t\t\t\t: %g\n'],...
        p.timeWindow, numel(allCorLen), numel(finiteCorLen))];
        
    for j=1:numel(corLenValues)
        channelLog{i} = [channelLog{i} ...
           sprintf('Total number of points tracked with a template size of  %d\t: %g\n',...
        corLenValues(j),corLenCount(j))];  
    end
    clear allCorLen 
end

% Close waitbar
if ishandle(wtBar), close(wtBar); end
disp('Finished tracking flow!')

% Create process report
procLog=[sprintf('Flow tracking summary\n\n') channelLog{:}];
disp(procLog);
fid=fopen([p.OutputDirectory filesep 'FlowTrackingSummary.txt'],'wt');
fprintf(fid,'%s',procLog);
fclose(fid);