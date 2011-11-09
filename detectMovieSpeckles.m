function detectMovieSpeckles(movieData,varargin)
% detectMovieSpeckle detect speckles of a movie
%
% detectMovieSpeckles 
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
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous speckle detection processes                                                                     
iProc = movieData.getProcessIndex('SpeckleDetectionProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(SpeckleDetectionProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end
specDetProc = movieData.processes_{iProc};
%Parse input, store in parameter structure
p = parseProcessParams(specDetProc,paramsIn);

%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...','Name',specDetProc.getName());
else
    wtBar=-1;
end

% Reading various constants
imDirs  = movieData.getChannelPaths;
bitDepth = movieData.camBitdepth_;
nFrames = movieData.nFrames_;
maxIntensity =(2^bitDepth-1);

%Find the required segmentation process
iMaskProc =movieData.getProcessIndex('MaskRefinementProcess',1,1); 
if isempty(iMaskProc)
    error(['Mask refinement have not yet been performed '...
    'on this movie! Please run first!!']);
end

segProc = movieData.processes_{iMaskProc};
if ~all(segProc.checkChannelOutput(p.MaskChannelIndex))
    error(['Mask refinement has not been run for all selected channels! '...
        'Please apply segmentation before running correlation!']);
end
       

% Create mask directory if several masks need to be merged
if length(p.MaskChannelIndex) >1
    %Get the indices of any previous mask intersection process
    iMaskIntProc = movieData.getProcessIndex('MaskIntersectionProcess',1,0);
    
    %If the process doesn't exist, create it
    if isempty(iMaskIntProc)
        iMaskIntProc = numel(movieData.processes_)+1;
        movieData.addProcess(MaskIntersectionProcess(movieData,p.OutputDirectory));
    end
    maskIntProc = movieData.processes_{iMaskIntProc};
    
    %Set up the parameters for mask transformation
    maskIntParams.ChannelIndex = p.MaskChannelIndex;
    maskIntParams.SegProcessIndex = iMaskProc;
    
    parseProcessParams(maskIntProc,maskIntParams);
    maskIntProc.run;
    
    % Get mask directory and names
    maskDir = maskIntProc.outFilePaths_{1};
    maskNames = maskIntProc.getOutMaskFileNames(1);
else
    maskDir = segProc.outFilePaths_{p.MaskChannelIndex};
    maskNames = segProc.getOutMaskFileNames(p.MaskChannelIndex);
end

iNDProc =movieData.getProcessIndex('NoiseEstimationProcess',1,1);   
if ~isempty(iNDProc)
    noiseProc=movieData.processes_{iNDProc};
    if ~noiseProc.checkChannelOutput()
        error(['There is no noise estimation output !' ...
            'Please apply noise estimation before running speckle detection!'])
    end
    
    % Read noise parameters from noise estimation process if input
    [p.I0 p.sDN p.GaussRatio] = noiseProc.loadChannelOutput('output',{'I0','sDN','GaussRatio'});

end

% Set up the input directories
inFilePaths = cell(2,numel(movieData.channels_));
for j = p.ChannelIndex
    inFilePaths{1,j} = imDirs{j};
    inFilePaths{2,j} = maskDir;
    if ~isempty(iNDProc)
        inFilePaths{3,j} = movieData.processes_{iNDProc}.outFilePaths_{j};
    end
end
specDetProc.setInFilePaths(inFilePaths);
    
% Set up the output directories
outputDir = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex;    
    %Create string for current directory
    outputDir{i} = [p.OutputDirectory filesep 'speckles_for_channel_' num2str(i)];
    mkClrDir(outputDir{i});
end
specDetProc.setOutFilePaths(outputDir);

%% --------------- Speckle detection ---------------%%% 

disp('Starting detecting speckles...')
%Format string for zero-padding file names
fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];
numStr = @(frame) num2str(frame,fString);

% Anonymous functions for reading input/output
outFile=@(chan,frame) [outputDir{chan} filesep 'cands_' numStr(frame) '.mat'];
inMask=@(frame) [maskDir filesep maskNames{1}{frame}];

logMsg = @(chan) ['Please wait, detecting speckles for channel ' num2str(chan)];
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t)) 's'];
tic;
nChan = length(p.ChannelIndex);
nTot = nChan*nFrames;

channelLog=cell(1,numel(p.ChannelIndex));

for i = 1:numel(p.ChannelIndex)
    iChan = p.ChannelIndex(i);
    % Log display
    channelLog{i} = sprintf('Channel %g: %s\n',iChan,imDirs{iChan});
    disp(logMsg(iChan))
    disp(imDirs{iChan});
    disp('Results will be saved under:')
    disp(outputDir{iChan});
    if ishandle(wtBar), waitbar(0,wtBar,logMsg(iChan)); end

    % Construct noiseParameters vector from noise model and alpha value
    k = fzero(@(x)diff(normcdf([-Inf,x]))-1+p.alpha,1);
    noiseParam = [k/p.GaussRatio p.sDN 0 p.I0];
    
    allcands=[];
    for j= 1:nFrames
        % Load the current image, scale it and apply Gaussian filter
        currImage = movieData.channels_(iChan).loadImage(j)/maxIntensity; 
        if p.filterSigma(iChan)>0
            currImage = filterGauss2D(currImage,p.filterSigma(iChan));
        end
        
        % Mask the filtered image
        currImage= currImage.*logical(imread(inMask(j)));

        % Statistically test the local maxima to extract (significant) speckles
        [cands locMax] = detectSpeckles(currImage,noiseParam,p.paramSpeckles,p.filterSigma(iChan));  %#ok<NASGU>
        
        % Save results
        save(outFile(iChan,j), 'cands','locMax');
        allcands = cat(1,allcands,cands);
        
        % Update the 
        if mod(j,5)==1 && ishandle(wtBar)
            tj=toc;
            nj = (i-1)*nFrames+ j;
            waitbar(nj/nTot,wtBar,sprintf([logMsg(iChan) timeMsg(tj*nTot/nj-tj)]));
        end
    end
    
    % Create channel log fot output
    channelLog{i} = [channelLog{i} ...
        sprintf(['Total number of detected speckles in the movie\t\t\t: %g\n'...
        'Average number of detected speckles per frame\t\t\t: %g\n'...
        'Average number of statistically significant speckles per frame\t: %g\n'],...
        numel(allcands),numel(allcands)/nFrames,sum([allcands.status])/nFrames)];
    
    order = {'primary','secondary','tertiary','quaternary'};
    for j=1:max([allcands.speckleType])
        channelLog{i} = [channelLog{i} ...
           sprintf(['Average number of ' order{j} ' speckles per frame \t\t\t: %g\n'],...
        sum([allcands.speckleType]==j)/nFrames)];  
    end
    clear allcands
end
% Close waitbar
if ishandle(wtBar), close(wtBar); end

disp('Finished detecting speckles!')

% Create process report
procLog=[sprintf('Speckle detection summary\n\n') channelLog{:}];
disp(procLog);
fid=fopen([p.OutputDirectory filesep 'SpeckleDetectionSummary.txt'],'w+');
fprintf(fid,procLog);
fclose(fid);