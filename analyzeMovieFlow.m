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
imDirs  = movieData.getChannelPaths;
imageFileNames = movieData.getImageFileNames;
nFrames = movieData.nFrames_;

% Test the presence and output validity of the speckle tracking process
iSpecTrack =movieData.getProcessIndex('SpeckleTrackingProcess',1,1);     
if isempty(iSpecTrack)
    error(['Speckle tracking has not yet been performed'...
    'on this movie! Please run first!!']);
end        
%Check that there is a valid output
specTrackProc = movieData.processes_{iSpecTrack};
if ~specTrackProc.checkChannelOutput(p.ChannelIndex)
    error(['Each channel must have tracks!' ...
        'Please apply speckle tracking to all needed channels before '...
        'running flow analysis!'])
end
    
iSegProc =movieData.getProcessIndex('MaskRefinementProcess',1,1);     
if isempty(iSegProc)
    error(['Segmentation has not yet been performed'...
    'on this movie! Please run first!!']);
end        
%Check that there is a valid output
segProc = movieData.processes_{iSegProc};
if ~segProc.checkChannelOutput(p.ChannelIndex)
    error(['Each channel must have masks !' ...
        'Please apply speckle tracking to all needed channels before'...
        'running speckle tracking!'])
end

% Set up the input directories
inFilePaths = cell(2,numel(movieData.channels_));
for j = p.ChannelIndex
    inFilePaths{1,j} = specTrackProc.outFilePaths_{1,j};
    inFilePaths{2,j} = segProc.outFilePaths_{1,j};
end
flowAnProc.setInFilePaths(inFilePaths);
    
% Set up the output directories
outputDir=cell(1,numel(movieData.channels_));
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
inImage=@(chan,frame) [imDirs{chan} filesep imageFileNames{chan}{frame}];
logMsg = @(chan) ['Please wait, analyzing flow for channel ' num2str(chan)];
outFile=@(chan,frame) [outputDir{chan} filesep 'flowMaps_' numStr(frame) '.mat'];

for iChan = p.ChannelIndex
    % Log display
    disp(logMsg(iChan))
    
    if ishandle(wtBar), waitbar(0,wtBar,'Loading masks and images...'); end
    maskNames = segProc.getOutMaskFileNames(iChan);
    inMask=@(frame) [flowAnProc.inFilePaths_{2,iChan}...
        filesep maskNames{1}{frame}];
    mask=true([movieData.imSize_ nFrames]);
    stack=zeros([movieData.imSize_ nFrames]);
    for j = 1:nFrames
        mask(:,:,j) = logical(imread(inMask(j)));
        stack(:,:,j) = double(imread(inImage(iChan,j)));
    end
    
   
    % Load candidates and generate Nx3 matrices with position and intensity
    % Replace fsmTrackFillSpeckleList
    if ishandle(wtBar), waitbar(0,wtBar,'Loading tracks...'); end
    M = specTrackProc.loadChannelOutput(iChan,'output','M'); 
    replicateFrames = @(x) [repmat(x(1),1,fix(p.timeWindow/2)) x ...
        repmat(x(end),1,fix(p.timeWindow/2)+1)];
    
    % Interpolate field
    if ishandle(wtBar), waitbar(.25,wtBar,'Interpolating flow...'); end
    [Mv,Md,Ms,E,S] = ...
        analyzeFlow(M,p.timeWindow,p.corrLength,...
        'interpolate',p.interpolate,'noise',p.noise,'error',p.error);
    Md=replicateFrames(Md); %#ok<NASGU>
    Ms=replicateFrames(Ms); %#ok<NASGU>
    E=replicateFrames(E);
    S=replicateFrames(S);
    
    
    % Speed maps creation
    if ishandle(wtBar), waitbar(.5,wtBar,'Generating speed maps...'); end
    % Interpolate raw vector on a grid
    G=framework(movieData.imSize_,[p.gridSize p.gridSize]);
    Mdgrid=arrayfun(@(i) vectorFieldAdaptInterp(Mv{i},G,p.corrLength,...
        [],'strain'),1:size(M,3),'UniformOutput',false);

    speedMap = createSpeedMaps(cat(3,Mdgrid{:}),p.timeWindow,movieData.timeInterval_,...
        movieData.pixelSize_,movieData.imSize_,mask);
    speedMap=replicateFrames(speedMap); %#ok<NASGU>
    
    if ishandle(wtBar), waitbar(.75,wtBar,'Generating error maps...'); end
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
end
% Close waitbar
if ishandle(wtBar), close(wtBar); end

disp('Finished analyzing flow!');
