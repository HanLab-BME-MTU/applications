function analyzeIntensityModesOO(movieData,varargin)
%analyzeIntensityModesOO is a wrapper for analyzeIntensityModes for the MovieData framework
%
%Khuloud Jaqaman, March 2015

%% Input

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn = ip.Results.paramsIn;

%Get the indices of any previous intensity modal analysis processes                                                                         
%If the process doesn't exist, create it
iProc = movieData.getProcessIndex('IntModalAnalysisProcess',1,0);
if isempty(iProc)
    iProc=numel(movieData.processes_)+1;
    movieData.addProcess(IntModalAnalysisProcess(movieData));
end
modProc = movieData.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(modProc,paramsIn);

%% --------------- Initialization ---------------%%

% Get detection process index
iLocProc = movieData.getProcessIndex('DetectionProcess',1,1);

assert(~isempty(iLocProc),['Particles have not been detected! '...
    'Please run detection prior to intensity modal analysis!']);
locProc=movieData.processes_{iLocProc};

assert(all(locProc.checkChannelOutput(p.ChannelIndex)),...
    ['Particles have not been detected! Please run detection before ' ...
    'running  intensity modal analysis!']);

% Set up the input directories (detection)
inFilePaths = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex
    inFilePaths{1,i} = locProc.outFilePaths_{1,i};
end
modProc.setInFilePaths(inFilePaths);

% Set up the output file
outFilePaths = cell(1,numel(movieData.channels_));
mkClrDir(p.OutputDirectory)
for i = p.ChannelIndex;
    outFilePaths{1,i} = [p.OutputDirectory filesep 'channel_' num2str(i) '.mat'];
end
modProc.setOutFilePaths(outFilePaths);

%% --------------- intensity modal analysis ---------------%%%

disp('Intensity modal analysis ...');

for i = p.ChannelIndex
    
    %read input
    movieInfo = locProc.loadChannelOutput(i);
    
    %define figure name and file for saving
    movieName = [movieData.movieDataFileName_(1:end-4) '_channel_' num2str(i)];
    figFileName = [p.OutputDirectory filesep 'channel_' num2str(i) '_IntensityModalAnalysis.fig'];
    
    %call function to do modal analysis
    [fitResults,numModes,allModeMean,allModeStd,allModeFrac,...
    numFeatures] = analyzeIntensityModes(movieInfo,movieName,...
    p.startFrame,p.endFrame,p.alpha,p.variableMean,p.variableStd,...
    p.numModeMinMax,p.plotResults,p.logData,p.modeParamIn,p.ampOrInt,figFileName); %#ok<ASGLU>

    %save output
    save(outFilePaths{1,i} ,'fitResults','numModes','allModeMean',...
        'allModeStd','allModeFrac','numFeatures');
    
end

disp('Done');

