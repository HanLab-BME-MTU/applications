function analyzeMovieSpeckles(movieData,varargin)
%analyzeMovieSpeckles analyzes the tracked speckles and create kinetic maps
%
% SYNOPSIS analyzeMovieSpeckles(movieData,paramsIn)
%
% INPUT   
%   movieData - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%
% OUTPUT   

% Sebastien Besson, July 2011

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous speckle detection processes                                                                     
iProc = movieData.getProcessIndex('KineticAnalysisProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(KineticAnalysisProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end

kinProc = movieData.processes_{iProc};
%Parse input, store in parameter structure
p = parseProcessParams(kinProc,paramsIn);

%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...','Name',kinProc.getName());
    wtBarArgs={'waitbar',wtBar};
else
    wtBar=-1;
    wtBarArgs={};
end

% Reading various constants
imDirs  = movieData.getChannelPaths;
imageFileNames = movieData.getImageFileNames;
nFrames = movieData.nFrames_;
bitDepth = movieData.camBitdepth_;
maxIntensity =(2^bitDepth-1);

% Test the presence and output validity of the speckle tracking process
iSpecTrack =movieData.getProcessIndex('SpeckleTrackingProcess',1,1);     
if isempty(iSpecTrack)
    error(['Speckle tracking has not yet been performed'...
    'on this movie! Please run first!!']);
end        
%Check that there is a valid output
specTrackProc = movieData.processes_{iSpecTrack};
if ~specTrackProc.checkChannelOutput(p.ChannelIndex)
    error(['Each channel must have tracks! ' ...
        'Please apply speckle tracking to all needed channels before '...
        'running kinetic analysis!'])
end


iSpecDet =movieData.getProcessIndex('SpeckleDetectionProcess',1,1);     
specDetProc = movieData.processes_{iSpecDet};
p.alpha = specDetProc.funParams_.alpha;

% Load optional noise estimation process output
iNDProc =movieData.getProcessIndex('NoiseEstimationProcess',1,1);   
if ~isempty(iNDProc)
    noiseProc=movieData.processes_{iNDProc};
    if ~noiseProc.checkChannelOutput(p.ChannelIndex)
        error(['There is no noise estimation output !' ...
            'Please apply noise estimation before running speckle detection!'])
    end

else
    % Read noise parameters from speckle detection process
    p.I0=specDetProc.funParams_.I0;
    p.sDN=specDetProc.funParams_.sDN;
    p.GaussRatio=specDetProc.funParams_.GaussRatio;
end


% Set up the input directories
nChan=numel(movieData.channels_);
inFilePaths = cell(2,nChan);
for j = p.ChannelIndex
    inFilePaths{1,j} = specDetProc.outFilePaths_{1,j};
    inFilePaths{2,j} = specTrackProc.outFilePaths_{1,j};
    if ~isempty(iNDProc)
        inFilePaths{3,j} = movieData.processes_{iNDProc}.outFilePaths_{1,j};
    end
end
kinProc.setInFilePaths(inFilePaths);
    
% Set up the output directories
outputDir=cell(1,nChan);
for i = p.ChannelIndex;    
    %Create string for current directory
    outputDir{i} = [p.OutputDirectory filesep 'channel_' num2str(i)];
    mkClrDir(outputDir{i});
end
kinProc.setOutFilePaths(outputDir);

%% --------------- Kinetic analysi ---------------%%% 

disp('Starting analyzing speckles...')
%Format string for zero-padding file names
fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];
numStr = @(frame) num2str(frame,fString);
% Anonymous functions for reading input/output
inImage=@(chan,frame) [imDirs{chan} filesep imageFileNames{chan}{frame}];
logMsg = @(chan) ['Please wait, analyzing speckles for channel ' num2str(chan)];
outFile=@(chan,frame) [outputDir{chan} filesep 'kineticMaps_' numStr(frame) '.mat'];

kineticLimits=cell(1,nChan);
channelLog=cell(1,numel(p.ChannelIndex));

for i = 1:numel(p.ChannelIndex);
    iChan = p.ChannelIndex(i);
    % Log display
    channelLog{i} = sprintf('Channel %g: %s\n',iChan,imDirs{iChan});
    disp(logMsg(iChan))
    if ishandle(wtBar), set(wtBar,'Name',[kinProc.getName() ' - Channel ' num2str(iChan)]); end
    
    % Load images
    fprintf(1,'Loading images...\n'); 
    filteredImages = cell(1,nFrames);
    for j= 1:nFrames        
        %Load the current image and apply Gaussian filter
        filteredImages{j} = double(imread(inImage(iChan,j)))/maxIntensity;
        filterSigma = specDetProc.funParams_.filterSigma(iChan);
        if filterSigma>0
            filteredImages{j} = filterGauss2D(filteredImages{j},filterSigma);
        end
    end

    % Load candidates and generate Nx3 matrices with position and intensity
    % Replace fsmTrackFillSpeckleList
    fprintf(1,'Loading speckles...\n');    
    cands = specDetProc.loadChannelOutput(iChan);
    
    % Load candidates and generate Nx3 matrices with position and intensity
    % Replace fsmTrackFillSpeckleList
    fprintf(1,'Loading tracks...\n');    
    [MPM,gapList] = specTrackProc.loadChannelOutput(iChan,'output',{'MPM','gapList'});

    % Read noise parameters from noise estimation process if input
    if ~isempty(iNDProc)
         [p.I0 p.sDN p.GaussRatio] = noiseProc.loadChannelOutput(iChan,...
             'output',{'I0','sDN','GaussRatio'});
    end

    % Create noise parameter vector
    k = fzero(@(x)diff(normcdf([-Inf,x]))-1+p.alpha,1);
    noiseParam = [k/p.GaussRatio p.sDN 0 p.I0];
    threshold = specTrackProc.funParams_.threshold;
    
    speckleArray=buildSpeckleArray(MPM,noiseParam,threshold,gapList,cands,...
        filteredImages,wtBarArgs{:});
    [speckleArray,log]=classifyKineticEvents(speckleArray,p.bleachRed,k,wtBarArgs{:});
    channelLog{i} = [channelLog{i} log];
    
    % Get activity, birth and death events index
    scoreIndx=[speckleArray.activity]~=0;
    birthEvents=speckleArray.status(scoreIndx)=='b';
    deathEvents=speckleArray.status(scoreIndx)=='d';

    % Construct double score matrix
    score=[double(speckleArray.timepoint(scoreIndx)) ...
        double(speckleArray.spPos(scoreIndx,:)) ...
        speckleArray.score(scoreIndx)];
    % Accounts for births and deaths in timepoints data
    score(birthEvents,1)=score(birthEvents,1)+1;
    score(deathEvents,1)=score(deathEvents,1)-1;
    % Sort by timepoint
    score=sortrows(score,1:4);

    % Construct kinScore cell array
    kinScore = arrayfun(@(x)score(score(:,1)==x,:),1:movieData.nFrames_,...
        'UniformOutput',false);
    noKinScore = cellfun(@isempty,kinScore);
    kinScoreArray = [find(noKinScore)' zeros(sum(noKinScore),3)];
    kinScore(noKinScore)=num2cell(kinScoreArray,2)';
    
    % Create averaged kinetic maps and replicate maps
    replicateFrames = @(x) [repmat(x(1),1,fix(p.timeWindow/2)) x ...
        repmat(x(end),1,fix(p.timeWindow/2))];
    [polyMap,depolyMap,kinMap2C] = createKineticMaps(kinScore,...
        p.timeWindow,movieData.imSize_,p.sigma,wtBarArgs{:});
    polyMap=replicateFrames(polyMap); 
    depolyMap=replicateFrames(depolyMap);
    kinMap2C=replicateFrames(kinMap2C);
    
    % Save output
    disp('Results will be saved under:')
    disp(kinProc.outFilePaths_{1,iChan});
    for j=1:nFrames
        s.polyMap=polyMap{j};
        s.depolyMap=depolyMap{j};
        s.kinMap2C=kinMap2C{j};
        save(outFile(iChan,j),'-struct','s');
    end
    
    % Store kinetic map limits
    allPolyMaps=vertcat(polyMap{:});
    allDepolyMaps=vertcat(depolyMap{:});
    kineticLimits{iChan}=[min(allDepolyMaps(:)) max(allPolyMaps(:))];
end
kinProc.setKineticLimits(kineticLimits)

% Close waitbar
if ishandle(wtBar), close(wtBar); end

disp('Finished analyzing speckles!');

% Create process report
procLog=[sprintf('Kinetic analysis summary\n\n') channelLog{:}];
disp(procLog);
fid=fopen([p.OutputDirectory filesep 'KineticAnalysisSummary.txt'],'w+');
fprintf(fid,procLog);
fclose(fid);

