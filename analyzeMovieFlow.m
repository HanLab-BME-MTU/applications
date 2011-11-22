function analyzeMovieFlow(movieData,varargin)
% analyzeMovieFlow reports statistics on vector fields
%
% SYNOPSIS analyzeMovieFlow(movieData,paramsIn)
%
% INPUT   
%   movieData - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%
% OUTPUT   

% Sebastien Besson, June 2011

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous speckle detection processes                                                                     
iProc = movieData.getProcessIndex('FlowAnalysisProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(FlowAnalysisProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end

flowAnProc = movieData.processes_{iProc};
%Parse input, store in parameter structure
p = parseProcessParams(flowAnProc,paramsIn);

%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...','Name',flowAnProc.getName());
end

% Reading various constants
nFrames = movieData.nFrames_;

% Test the presence and output validity of the speckle tracking process
iSpecDetProc =movieData.getProcessIndex('SpeckleDetectionProcess',1,1);     
if isempty(iSpecDetProc)
    error(['Speckle detection has not yet been performed'...
    'on this movie! Please run first!!']);
end        
%Check that there is a valid output
specDetProc = movieData.processes_{iSpecDetProc};
if ~specDetProc.checkChannelOutput(p.ChannelIndex)
    error(['Each channel must have speckles! ' ...
        'Please apply speckle detection to all needed channels before '...
        'running flow analysis!'])
end
p.MaskChannelIndex = specDetProc.funParams_.MaskChannelIndex;
    

% Create mask directory if several masks need to be merged
if length(p.MaskChannelIndex) >1
    %Get the indices of any previous mask intersection process
    iMaskProc = movieData.getProcessIndex('MaskIntersectionProcess',1,0);
else  
    iMaskProc =movieData.getProcessIndex('MaskRefinementProcess',1,1);
end

if isempty(iMaskProc)
    error(['Mask refinement has not yet been performed '...
        'on this movie! Please run first!!']);
end
%Check that there is a valid output
maskProc = movieData.processes_{iMaskProc};
if ~maskProc.checkChannelOutput(p.ChannelIndex)
    error(['Each channel must have masks !' ...
        'Please apply mask refinement to all needed channels before'...
        'running speckle tracking!'])
end


% Test the presence and output validity of the speckle tracking process
iFlowProc =movieData.getProcessIndex(p.FlowProcess,1,1);     
if isempty(iFlowProc)
    error([eval([p.FlowProcess '.getName']) ' has not yet been performed'...
    'on this movie! Please run first!!']);
end        

%Check that there is a valid output
flowProc = movieData.processes_{iFlowProc};
if ~flowProc.checkChannelOutput(p.ChannelIndex)
    error(['Each channel must have flow! Please apply '...
        eval([p.FlowProcess '.getName']) ' to all needed channels before '...
        'running flow analysis!'])
end
    

% Set up the input directories
nChan=numel(movieData.channels_);
inFilePaths = cell(2,nChan);
for j = p.ChannelIndex
    inFilePaths{1,j} = flowProc.outFilePaths_{1,j};
    inFilePaths{2,j} = maskProc.outFilePaths_{1,j};
end
flowAnProc.setInFilePaths(inFilePaths);
    
% Set up the output directories
outputDir=cell(1,nChan);
for i = p.ChannelIndex;    
    %Create string for current directory
    outputDir{i} = [p.OutputDirectory filesep 'channel_' num2str(i)];
    mkClrDir(outputDir{i});
end
flowAnProc.setOutFilePaths(outputDir);

%% --------------- Kinetic analysi ---------------%%% 

disp('Starting analyzing flow...')
%Format string for zero-padding file names
fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];
numStr = @(frame) num2str(frame,fString);
% Anonymous functions for reading input/output
logMsg = @(chan) ['Please wait, analyzing flow for channel ' num2str(chan)];
outFile=@(chan,frame) [outputDir{chan} filesep 'flowMaps_' numStr(frame) '.mat'];

speedMapLimits=cell(1,nChan);
flowLimits=cell(1,nChan);
for iChan = p.ChannelIndex
    % Log display
    disp(logMsg(iChan))
    
    if ishandle(wtBar), waitbar(0,wtBar,'Loading masks and images...'); end
    maskNames = maskProc.getOutMaskFileNames(iChan);
    inMask=@(frame) [flowAnProc.inFilePaths_{2,iChan}...
        filesep maskNames{1}{frame}];
    mask=true([movieData.imSize_ nFrames]);
    stack=zeros([movieData.imSize_ nFrames]);
    for j = 1:nFrames
        mask(:,:,j) = logical(imread(inMask(j)));
        stack(:,:,j) = movieData.channels_(iChan).loadImage(j);
    end
    
   
    % Load candidates and generate Nx3 matrices with position and intensity
    % Replace fsmTrackFillSpeckleList
    if ishandle(wtBar), waitbar(0,wtBar,['Loading flow for channel ' num2str(iChan)']); end
    
    switch p.FlowProcess
        case 'SpeckleTrackingProcess'
            M = flowProc.loadChannelOutput(iChan,'output','M');
            flow=arrayfun(@(i)M(M(:,1,i)~=0 & M(:,3,i)~=0,:,i),1:size(M,3),'Unif',false);
        case 'FlowTrackingProcess'
            flow = flowProc.loadChannelOutput(iChan,'output','flow');
            flow=flow(1:end-1);
    end
    
    % Interpolate field
    if ishandle(wtBar), waitbar(.25,wtBar,['Interpolating flow for channel ' num2str(iChan)']); end
    [Md,Ms,E,S] = ...
        analyzeFlow(flow,p.timeWindow,p.corrLength,'noise',p.noise,'error',p.error);
    
    
    % Repliacte frames
    replicateFrames = @(x) [repmat(x(1),1,fix(p.timeWindow/2)) x ...
        repmat(x(end),1,fix(p.timeWindow/2)+1)];
    Md=replicateFrames(Md); 
    Ms=replicateFrames(Ms); %#ok<NASGU>
    E=replicateFrames(E);
    S=replicateFrames(S);
    
    
    % Speed maps creation
    if ishandle(wtBar), waitbar(.5,wtBar,['Generating speed maps for channel ' num2str(iChan)']); end
    % Interpolate raw vector on a grid
    G=framework(movieData.imSize_,[p.gridSize p.gridSize]);
    Mdgrid=arrayfun(@(i) vectorFieldAdaptInterp(flow{i},G,p.corrLength,...
        [],'strain'),1:size(M,3),'UniformOutput',false);

    speedMap = createSpeedMaps(cat(3,Mdgrid{:}),p.timeWindow,movieData.timeInterval_,...
        movieData.pixelSize_,movieData.imSize_,mask);
    speedMap=replicateFrames(speedMap);
    
    if ishandle(wtBar), waitbar(.75,wtBar,['Generating error maps for channel ' num2str(iChan)']); end
    [img3C_map img3C_SNR]=createErrorMaps(stack,E,S); %#ok<ASGLU,NASGU>
    
    % Fill output structure for each frame and save it
    disp('Results will be saved under:')
    disp(flowAnProc.outFilePaths_{1,iChan});
    output={'speedMap','Md','Ms','E','S','img3C_map','img3C_SNR'};
    for j=1:nFrames
        for k=1:numel(output)
            s.(output{k})=eval([output{k} '{' num2str(j) '}']);
        end
        save(outFile(iChan,j),'-struct','s');
    end
    
    % Store speed maps and flow limits
    allMaps = vertcat(speedMap{:});
    speedMapLimits{iChan}=[min(allMaps(:)) max(allMaps(:))];
    
    allFlow = vertcat(Md{:});
    flowLimits{iChan}=[min(allMaps(:)) max(allMaps(:))];
    flowMagnitude = (diff(allFlow(:,[1 3]),1,2).^2+diff(allFlow(:,[2 4]),1,2).^2).^.5;
    flowLimits{iChan}=[min(flowMagnitude(:)) max(flowMagnitude(:))];
    
end
flowAnProc.setSpeedMapLimits(speedMapLimits)
flowAnProc.setFlowLimits(flowLimits);

% Close waitbar
if ishandle(wtBar), close(wtBar); end

disp('Finished analyzing flow!');
