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
wtBar = -1;
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
channelName = @(x)movieData.getChannelPaths{x}(max(regexp(movieData.getChannelPaths{x},filesep))+1:end);   
for i = p.ChannelIndex;    
    %Create string for current directory
    outputDir{i} = fullfile(p.OutputDirectory,channelName(i));
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
protSpeedMapLimits=cell(1,nChan);
flowLimits=cell(1,nChan);
channelLog=cell(1,numel(p.ChannelIndex));

for i=1:numel(p.ChannelIndex)
    iChan = p.ChannelIndex(i);
    % Log display
    channelLog{i} = sprintf('Channel %g: %s\n',iChan,inFilePaths{1,iChan});
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
            
            % Remove NaNs 
            for j=find(~cellfun(@isempty,flow))
                flow{j}=flow{j}(~isnan(flow{j}(:,3)),:);
            end           
    end
    
    % Drift correction
    if p.driftCorrection
        velField = cellfun(@(x) [x(:,4)-x(:,2),x(:,3)-x(:,1)],flow,'unif',false);
        meanVel = cellfun(@(x) mean(x),velField,'unif',false);
        %Subtract meanVel from velField
        velFieldCorrected = cellfun(@(x,y) [x(:,1)-y(1), x(:,2)-y(2)],velField,meanVel,...
            'unif',false);
%         figure, quiver(flow{1}(:,2),flow{1}(:,1),...
%             flow{1}(:,4)-flow{1}(:,2),flow{1}(:,3)-flow{1}(:,1))
%         hold on
%         quiver(flow{1}(:,2),flow{1}(:,1),...
%             velFieldCorrected{1}(:,1),velFieldCorrected{1}(:,2))
        % Put it back to flow
        for ii=1:numel(flow)
            flow{ii}(:,3:4) = flow{ii}(:,1:2)+velFieldCorrected{ii}(:,2:-1:1); %#ok<NASGU>
        end
    end
    
    % Interpolate field
    if ishandle(wtBar), waitbar(.25,wtBar,['Interpolating flow for channel ' num2str(iChan)']); end
    [Md,Ms,E,S,stats] =  analyzeFlow(flow,p.timeWindow,p.corrLength,...
        'noise',p.noise,'error',p.error);

    % Filter vectors for outliers
    outlierThreshold = 10;
    for j=1:nFrames
        mdMag = ((Md{j}(:,3)-Md{j}(:,1)).^2+(Md{j}(:,4)-Md{j}(:,2)).^2).^0.5;
        outlierIndex = detectOutliers(mdMag,outlierThreshold);
        Md{j}(outlierIndex,:)=[];
        Ms{j}(outlierIndex,:)=[];
        E{j}(outlierIndex,:)=[];
        S{j}(outlierIndex,:)=[];
    end
    
    % Create speed maps
    if ishandle(wtBar), waitbar(.5,wtBar,['Generating speed maps for channel ' num2str(iChan)']); end
    speedMap = createSpeedMaps(flow,p.timeWindow,p.corrLength,movieData.timeInterval_,...
        movieData.pixelSize_,movieData.imSize_,p.gridSize,mask);
    
    % Create error mapss
    if ishandle(wtBar), waitbar(.6,wtBar,['Generating error maps for channel ' num2str(iChan)']); end
    [img3C_map img3C_SNR]=createErrorMaps(stack,E,S); %#ok<ASGLU,NASGU>

    % Create protrusive flow speed (with sign)
    if ishandle(wtBar), waitbar(.75,wtBar,['Generating protrusive speed maps for channel ' num2str(iChan)']); end
    protSpeedMap = createProtrusiveSpeedMaps(Md,mask,movieData.pixelSize_,movieData.timeInterval_,p.gridSize); 
        
    % Fill output structure for each frame and save it
    disp('Results will be saved under:')
    disp(flowAnProc.outFilePaths_{1,iChan});
    for j=1:nFrames
        s.Md=Md{j};
        s.Ms=Ms{j};
        s.E=E{j};
        s.S=S{j};
        s.speedMap=speedMap{j};
        s.img3C_map=img3C_map{j};
        s.img3C_SNR=img3C_SNR{j};
        s.protSpeedMap=protSpeedMap{j};
        
        save(outFile(iChan,j),'-struct','s');
    end
    
    % Store speed maps and flow limits
    allMaps = vertcat(speedMap{:});
    speedMapLimits{iChan}=[min(allMaps(:)) max(allMaps(:))];
    
    % Store speed maps and flow limits
    allMapsProt = vertcat(protSpeedMap{:});
    protSpeedMapLimits{iChan}=[min(allMapsProt(:)) max(allMapsProt(:))];

    allFlow = vertcat(Md{:});
    
    flowLimits{iChan}=[min(allMaps(:)) max(allMaps(:))];
    flowMagnitude = (diff(allFlow(:,[1 3]),1,2).^2+diff(allFlow(:,[2 4]),1,2).^2).^.5;
    flowLimits{iChan}=[min(flowMagnitude(:)) max(flowMagnitude(:))];
    
    % Create channel log fot output
    lv=vertcat(stats.lv{:})*(60/movieData.timeInterval_)*movieData.pixelSize_;
    ld=vertcat(stats.ld{:})*(60/movieData.timeInterval_)*movieData.pixelSize_;
    ls=vertcat(stats.ls{:})*(60/movieData.timeInterval_)*movieData.pixelSize_;
    snr=vertcat(stats.snr{:});
    
    channelLog{i} = [channelLog{i} ...
        sprintf('Number of RAW vectors            : %d\n',numel(lv))...
        sprintf('Mean RAW vector length           : %2.4f nm/min +/- %2.4f nm/min (+/- %2.2f%%)\n',mean(lv),std(lv),100*std(lv)/mean(lv))...
        sprintf('Median RAW vector length         : %2.4f nm/min \n',median(lv))...
        sprintf('Mean INTERPOLATED vector length  : %2.4f nm/min +/- %2.4f nm/min (+/- %2.2f%%)\n',mean(ld),std(ld),100*std(ld)/mean(ld))...
        sprintf('Median INTERPOLATED vector length: %2.4f nm/min\n',median(ld))...
        sprintf('Mean NOISE vector length         : %2.4f nm/min +/- %2.4f nm/min (+/- %2.2f%%)\n',mean(ls),std(ls),100*std(ls)/mean(ls))...
        sprintf('Median NOISE vector length       : %2.4f nm/min\n',median(ls))...
        sprintf('Mean / median SNR                : %2.4f +/- %2.4f (+/- %2.2f%%) / %2.4f\n',mean(snr),std(snr),100*std(snr)/mean(snr),median(snr))];
    
end
flowAnProc.setSpeedMapLimits(speedMapLimits)
flowAnProc.setProtSpeedMapLimits(protSpeedMapLimits)
flowAnProc.setFlowLimits(flowLimits);

% Close waitbar
if ishandle(wtBar), close(wtBar); end

disp('Finished analyzing flow!');

% Create process report
procLog=[sprintf('Flow analysis summary\n\n') channelLog{:}];
disp(procLog);
fid=fopen([p.OutputDirectory filesep 'FlowAnalysisSummary.txt'],'wt');
fprintf(fid,'%s',procLog);
fclose(fid);
