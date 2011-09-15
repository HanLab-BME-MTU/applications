function correctMovieDisplacementField(movieData,varargin)
% correctMovieDisplacementField calculate the displacement field
%
% correctMovieDisplacementField 
%
% SYNOPSIS correctMovieDisplacementField(movieData,paramsIn)
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
iProc = movieData.getProcessIndex('DisplacementFieldCorrectionProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(DisplacementFieldCorrectionProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end
displFieldCorrProc = movieData.processes_{iProc};
%Parse input, store in parameter structure
p = parseProcessParams(displFieldCorrProc,paramsIn);

%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...','Name',displFieldProc.getName());
end

% Reading various constants
bitDepth = movieData.camBitdepth_;
nFrames = movieData.nFrames_;
maxIntensity =(2^bitDepth-1);

% Check displacement field process
iDisplFieldCalcProc =movieData.getProcessIndex('DisplacementFieldCalculationProcess',1,1);     
if isempty(iDisplFieldCalcProc)
    error(['Displacement field calculation has not been run! '...
        'Please run displacement field calculation prior to force field calculation!'])   
end

displFieldCalcProc=movieData.processes_{iDisplFieldProc};
if ~displFieldCalcProc.checkChannelOutput
    error(['The channel must have a displacement field ! ' ...
        'Please calculate displacement field to all needed channels before '...
        'running force field calculation!'])
end

inFilePaths{1} = displFieldCalcProc.outFilePaths_{1};
displFieldCorrProc.setInFilePaths(inFilePaths);

% Set up the output directories
outputFile = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex;    
    %Create string for current directory
    outputFile{i} = [p.OutputDirectory filesep 'displField.mat'];
    mkClrDir(p.OutputDirectory);
end
displFieldCorrProc.setOutFilePaths(outputFile);

%% --------------- Displacement field calculation ---------------%%% 

disp('Starting correctiong displacement field...')
% Anonymous functions for reading input/output
displField=displFieldCalcProc.loadChannelOutput;


disp('Detecting and filtering vector field outliers...')
logMsg = 'Please wait, detecting and filtering vector field outliers';
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t/60)) 'min'];
tic;

% Perform vector field outlier detection
if feature('ShowFigureWindows'), waitbar(0,wtBar,sprintf(logMsg)); end
for j= 1:nFrames

    % Outlier detection
    dispMat = [displField(j).pos displField(j).vec];
    if ~isempty(p.outlierThreshold)
        outlierIndex = detectVectorFieldOutliers(dispMat,p.outlierThreshold,1);
        displField(j).pos(outlierIndex,:)=[];
        displField(j).vec(outlierIndex,:)=[];
    end
    
    % Update the waitbar
    if mod(j,5)==1 && feature('ShowFigureWindows')
        tj=toc;
        waitbar(j/nFrames,wtBar,sprintf([logMsg timeMsg(tj*(nFrames-j)/j)]));
    end
end
% Find rotational registration
% if p.doRotReg
%    displField=perfRotReg(displField,1);
% end

save([p.OutputDirectory filesep 'displField.mat'],'displField');

% Close waitbar
if feature('ShowFigureWindows'), close(wtBar); end

disp('Finished correcting displacement field!')