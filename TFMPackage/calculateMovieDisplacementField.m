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
end

% Reading various constants
bitDepth = movieData.camBitdepth_;
nFrames = movieData.nFrames_;
maxIntensity =(2^bitDepth-1);

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
    imageFileNames = SDCProc.getInImageFileNames(p.ChannelIndex);
    s = load(SDCProc.outFilePaths_{3,p.ChannelIndex});
    residualT = s.T-round(s.T);
    refFrame = double(imread(SDCProc.outFilePaths_{2,p.ChannelIndex}));
else
    imDirs  = movieData.getChannelPaths(p.ChannelIndex);
    imageFileNames = movieData.getImageFileNames(p.ChannelIndex);
    refFrame = double(imread(p.referenceFramePath));
    residualT = zeros(nFrames,2);
end
inFilePaths{1,p.ChannelIndex} = imDirs{:};
displFieldProc.setInFilePaths(inFilePaths);
    
% Set up the output directories
outputFile = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex;    
    %Create string for current directory
    outputFile{i} = [p.OutputDirectory filesep 'displField.mat'];
    mkClrDir(p.OutputDirectory);
end
displFieldProc.setOutFilePaths(outputFile);

%% --------------- Displacement field calculation ---------------%%% 

disp('Starting calculating displacement field...')
% Anonymous functions for reading input/output
inImage=@(chan,frame) [imDirs{chan} filesep imageFileNames{chan}{frame}];

% Detect beads in reference frame and select only valid candidates
disp('Detecting beads in the reference frame...')
filteredRefFrame = filterGauss2D(refFrame/maxIntensity,...
    movieData.channels_(p.ChannelIndex(1)).psfSigma_);
k = fzero(@(x)diff(normcdf([-Inf,x]))-1+p.alpha,1);
noiseParam = [k/p.GaussRatio p.sDN 0 p.I0];
cands = detectSpeckles(filteredRefFrame,noiseParam,[1 0]);
M = vertcat(cands([cands.status]==1).Lmax);
beads = M(:,2:-1:1);

% For debugging purposes
% indx=beads(:,1)>600 & beads(:,1)<900&beads(:,2)>350&beads(:,2)<650;
% beads=beads(indx,:);

% Remove beads using the correlation length as the erosion distance
beadsMask = true(size(filteredRefFrame));
erosionDist=p.minCorLength+1;
beadsMask(erosionDist:end-erosionDist,erosionDist:end-erosionDist)=false;
indx=beadsMask(sub2ind(size(beadsMask),beads(:,2),beads(:,1)));
beads(indx,:)=[];

% Initialize displacement field structure
displField(nFrames)=struct('pos',[],'vec',[]);

disp('Calculating displacement field...')
logMsg = 'Please wait, calculating displacement field';
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t/60)) 'min'];
tic;

% Perform sub-pixel registration
for j= 1:1
    % Read image and perform correlation
    currImage = double(imread(inImage(p.ChannelIndex(1),j)));

    % Filter out beads which are in the background
    indx = currImage(sub2ind(size(currImage),beads(:,2),beads(:,1)))==0;
    localbeads = beads(~indx,:);
    
    [vx,vy] = trackStackFlow(cat(3,refFrame,currImage),...
        localbeads(:,1),localbeads(:,2),...
        p.minCorLength,p.minCorLength,'maxSpd',p.maxFlowSpeed);
    
    % Extract finite displacement and prepare displFiel
    validV = ~isnan(vx) & ~isnan(vy) & ~isinf(vx);
    displField(j).pos=[localbeads(validV,1) localbeads(validV,2)];
    displField(j).vec=[vx(validV)+residualT(j,1) vy(validV)+residualT(j,2)];
    
    % Update the waitbar
    if mod(j,5)==1 && feature('ShowFigureWindows')
        tj=toc;
        waitbar(j/nFrames,wtBar,sprintf([logMsg timeMsg(tj*(nFrames-j)/j)]));
    end
end

save([p.OutputDirectory filesep 'displField.mat'],'displField');

% Close waitbar
if feature('ShowFigureWindows'), close(wtBar); end

disp('Finished calculating displacement field!')