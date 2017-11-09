function [] = maskDetectedStructure(MD, varargin)
%maskDetectedStructure creates a mask that covers sub cellular structures detected by SubResolutionProcess
%
%SYNOPSIS function [] = maskDetectedStructure(MD, varargin)
%
%INPUT
%   MD : MovieData to which this process will be applied. MD should contain
%        SubResolutionProcess
%   funParams : parameter for this process
%       .channel         : The index of the channel of interest
%       .psfSigmaMult    : multiplicative factor used to determine the size
%                          of mask. radius = psfSigmaMult * pasfSigma
%       .outputDirectory : Directory where the mask (ouput) will be stored
%
%OUPUT
%
%NOTES
%   Mask array is in plot coordinate system [y, x]
%
%Tae H Kim, June 2015

%% Input
%Check input
ip = inputParser;
ip.CaseSensitive = false;
%standard process input
ip.addRequired('MD', @(x) isa(x,'MovieData'));
ip.addOptional('funParams',[],@isstruct);
ip.parse(MD, varargin{:});
%for easier calling
parameter = ip.Results.funParams;
if isempty(parameter)
    parameter = SubcellMaskProcess.getDefaultParams(MD.outputDirectory_);
end
channel = parameter.channel;
psfSigmaMult = parameter.psfSigmaMult;
outputDirectory = parameter.outputDirectory;
%Get the indices of any previous masking processes from this function                                                                              
iProc = MD.getProcessIndex('SubcellMaskProcess',1,0);
%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(MD.processes_)+1;
    MD.addProcess(SubcellMaskProcess(MD,MD.outputDirectory_));                                                                                                 
end
maskProcess = MD.processes_{iProc};
%% Initialization
% Check detection process first
detProcessIndex = MD.getProcessIndex('DetectionProcess',1,1);
if isempty(detProcessIndex)
    error(['Detection has not been run! '...
        'Please run detection prior to masking detected strucutures!'])
end
detProcess=MD.processes_{detProcessIndex};
if ~detProcess.checkChannelOutput(channel)
    error(['Missing detection output ! Please apply detection before ' ...
        'running masking detected strucutres!'])
end
%Obtains necessary data from detection process
movieFeatures = detProcess.loadChannelOutput(1);
psfSigma = detProcess.funParams_.detectionParam.psfSigma;
%easier variables to call
nFrames = length(movieFeatures);
radius = round(psfSigma .* psfSigmaMult);
xMax = MD.imSize_(2);
yMax = MD.imSize_(1);
% Set up the input directories (input images)
inFilePaths = cell(1,numel(MD.channels_));
for i = channel
    inFilePaths{1,i} = detProcess.outFilePaths_{1,i};
end
maskProcess.setInFilePaths(inFilePaths);
% Set up the output file
outFilePaths = cell(1,numel(MD.channels_));
for i = channel
    outFilename= ['Channel_' num2str(i) '_masking_result'];
    outFilePaths{1,i} = [outputDirectory filesep outFilename '.mat'];
end
mkClrDir(outputDirectory);
maskProcess.setOutFilePaths(outFilePaths);
% set values in maskProcess
maskProcess.imSize_ = [yMax, xMax];
maskProcess.nFrames_ = nFrames;
%% mask creation
mask = maskDetectedStructure_StandAlone(movieFeatures, nFrames, xMax, yMax, radius);
%% save
save(outFilePaths{1,channel},'mask');
MD.save;
end

