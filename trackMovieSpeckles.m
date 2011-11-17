function trackMovieSpeckles(movieData,varargin)
% trackMovieSpeckles tracks the speckles of a movie

%
% SYNOPSIS detectMovieSpeckles(movieData,paramsIn)
%
% INPUT   
%   movieData - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%
% OUTPUT   

% Sebastien Besson, May 2011 (last modified Sep 2011)

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addParamValue('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous speckle detection processes                                                                     
iProc = movieData.getProcessIndex('SpeckleTrackingProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(SpeckleTrackingProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end

specTrackProc = movieData.processes_{iProc};
%Parse input, store in parameter structure
p = parseProcessParams(specTrackProc,paramsIn);

%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows'),
    wtBar = waitbar(0,'Initializing...','Name',specTrackProc.getName());
    wtBarArgs={'waitbar',wtBar};
else
    wtBarArgs={};
end

% Reading various constants
nFrames = movieData.nFrames_;
nChan = numel(movieData.channels_);

% Test the presence and output validity of the speckle detection process
iSpecProc =movieData.getProcessIndex('SpeckleDetectionProcess',1,1);     
if isempty(iSpecProc)
    error(['Speckle detection has not yet been performed'...
    'on this movie! Please run first!!']);
end        
%Check that there is a valid output
specDetProc = movieData.processes_{iSpecProc};
if ~specDetProc.checkChannelOutput(p.ChannelIndex)
    error(['Each channel must have speckles !' ...
        'Please apply speckle detection to all needed channels before'...
        'running speckle tracking!'])
end
    
% Set up the input directories
inFilePaths = cell(1,nChan);
for j = p.ChannelIndex
    inFilePaths{1,j} = specDetProc.outFilePaths_{1,j};
end

% Check optional process Flow Tracking
iFTProc =movieData.getProcessIndex('FlowTrackingProcess',1,1);     
if ~isempty(iFTProc)
    flowTrackProc=movieData.processes_{iFTProc};
    if ~flowTrackProc.checkChannelOutput(p.ChannelIndex)
        error(['Each channel must have flow tracking output !' ...
            'Please apply flow tracking to all needed channels before'...
            'running speckle tracking!'])
    end
    for j = p.ChannelIndex
        inFilePaths{2,j} = flowTrackProc.outFilePaths_{1,j};
    end
    initCorLen = flowTrackProc.funParams_.maxCorLength;
else
    initCorLen = Inf;
end
specTrackProc.setInFilePaths(inFilePaths);
    
% Set up the output directories
outFilePaths_=cell(1,nChan);
mkClrDir(p.OutputDirectory)
for i = p.ChannelIndex;    
    %Create string for current directory
    outFilePaths_{1,i} = [p.OutputDirectory filesep 'tracks_for_channel_' num2str(i) '.mat'];
end
specTrackProc.setOutFilePaths(outFilePaths_);

%% --------------- Speckle detection ---------------%%% 

disp('Starting tracking speckles...')
% Anonymous functions for reading input/output
logMsg = @(chan) ['Please wait, tracking speckles for channel ' num2str(chan)];
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t)) 's'];
tic;
nTot = numel(p.ChannelIndex)*nFrames;

for iChan = p.ChannelIndex
    % Log display
    disp(logMsg(iChan))
    
    
    % Load candidates and generate Nx3 matrices with position and intensity
    % Replace fsmTrackFillSpeckleList
    if ishandle(wtBar), waitbar(0,wtBar,'Loading speckles...'); end
    fprintf(1,'Loading speckles...\n');    
    cands = specDetProc.loadChannelOutput(iChan);
    validCands = cellfun(@(x) x([x.status]==1),cands,'UniformOutput',false);
    speckles = cellfun(@(x) horzcat([vertcat(x.Lmax) vertcat(x.ILmax)]),...
        validCands,'UniformOutput',false);
    invalidCands = cellfun(@(x) x([x.status]==0),cands,'UniformOutput',false);
    clear validCands 
    
    % Initialize flow results to empty by default
    flow = cell(1,nFrames);
    if ~isempty(iFTProc)
        if ishandle(wtBar), waitbar(0,wtBar,'Loading flow field...'); end
        fprintf(1,'Loading flow field...');
        
        % Load flow field and correlation length for all frames
        flow = flowTrackProc.loadChannelOutput(iChan,'output','flow');
        
        % Interpolate non empty-flow
        for j=find(~cellfun(@isempty,flow))
            flow{j} = interpolateFlow(flow{j},p.corrLength,p.interpolationMethod);
        end
    end
    
    M=zeros(4,4,nFrames-1);
    vectorField = cell(1,nFrames-1);
    for j=1:nFrames-1;
        % Track speckles between frame j and j+1;
        if p.enhanced
            [matchM,vectors] = trackSpeckles(speckles{j},speckles{j+1},p.threshold,...
                'initM',flow{j},'initCorLen',initCorLen,...
                'enhanced',p.enhanced,'corrLength',p.corrLength);
            % Save interpolated vector field to disk for later use with gap closer
            vectorField{j}=vectors;
        else
            matchM = trackSpeckles(speckles{j},speckles{j+1},p.threshold,...
                'initM',flow{j},'initCorLen',initCorLen,...
                'enhanced',p.enhanced);
        end
        
        % If no flow tracked for frame j, use results of frame j
        if isempty(flow{j+1}) && ~isempty(iFTProc),
            % Extract vectors from M
            raw=matchM(matchM(:,1)~=0 & matchM(:,3)~=0,:);
            % Interpolate onto vector positions
            grid=raw(:,1:2);
            
            % Average returned M to be used to propagate I again
            interp_flow=vectorFieldAdaptInterp(raw,grid,p.corrLength,[],'strain');
            flow{j+1} = interp_flow;
        end
        M(1:size(matchM,1),1:size(matchM,2),j,1)=matchM;
        
        if mod(j,5)==1 && ishandle(wtBar)
            tj=toc;
            nj = sum(nFrames(1:i-1))+ j;
            waitbar(nj/nTot,wtBar,sprintf([logMsg(iChan) timeMsg(tj*nTot/nj-tj)]));
        end
    end

    % Correct dimensions
    cM=M>0; 
    [i,~,~]=find(cM); 
    M=M(1:max(i),:,:,:);
    
    % In case only two frames have been tracked, no gaps can be closed     
    if size(M,3)==1, MPM=M; end %#ok<NASGU>
    
    vector_arguments={};
    vector_save={};
    if p.enhanced,
        vector_arguments={'vectors',vectorField};
        vector_save={'vectorField'};
    end
    % Close gaps in M
    [M,gapList]=trackSpecklesGapCloser(M,p.threshold,invalidCands,...
        vector_arguments{:},wtBarArgs{:}); %#ok<NASGU>
    
    % Link gaps
    [MPM,M]=trackSpecklesLinker(M,wtBarArgs{:}); %#ok<NASGU,ASGLU>
    
    % Save output
    disp('Results will be saved as:')
    disp(specTrackProc.outFilePaths_{1,iChan});
    save(specTrackProc.outFilePaths_{1,iChan},'MPM','M','gapList','flow',...
        vector_save{:});
end
% Close waitbar
if ishandle(wtBar), close(wtBar); end

disp('Finished tracking speckles!');

function flow = interpolateFlow(flow,corLen, method)

% Find points needed to be interpolated
isnanIndx = isnan(flow(:,3)) | isnan(flow(:,4));
if isempty(isnanIndx), return; end

switch method
    case 'none' % Use vectors of magnitude 0
        flow(isnanIndx,3:4)=flow(isnanIndx,1:2);
    case 'nearest-neighbor'  % Use magnitude of nearest neighbor
        realFlow = flow(~isnanIndx,:);
        % Contruct Delaunay triangulation and find nearest neighbor
        % delaunayTri=DelaunayTri(realFlow(:,1),realFlow(:,2));
        % nnIdx = delaunayTri.nearestNeighbor(flow(isnanIndx,1),flow(isnanIndx,2));
        
        nnIdx = KDTreeClosestPoint(realFlow(:,1:2),flow(isnanIndx,1:2));
        flow(isnanIndx,3:4)=flow(isnanIndx,1:2)+...
            realFlow(nnIdx,3:4)-realFlow(nnIdx,1:2);
    case 'gaussian'  % Use flow tracking correlation length
        flow(isnanIndx,:) = vectorFieldInterp(flow(~isnanIndx,:),...
            flow(isnanIndx,1:2),corLen,[]);
end