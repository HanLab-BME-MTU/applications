function [] = trackPartitioning(MD_tracks, MD_struct, varargin)
%Determines the paritioning fraction of tracks into the mask
%
%SYNOPSIS function [] = trackPartitioning(MD_tracks, MD_struct, varargin)
%
%INPUT
%   MD_tracks : MovieData containing information about the tracks
%   MD_struct : MovieData containing information about the mask
%   funParams : Standard structure that contains the parameter for function
%               that handles the MovieData Processes. This function has no
%               parameters associated with it.
%
%OUTPUT
%
%Notes
%   Mask array is in plot coordinate system [y, x]
%
%Tae H Kim, June 2015
%% Input
%Check input
ip = inputParser;
ip.CaseSensitive = false;
%standard process input
ip.addRequired('MD_tracks', @(x) isa(x,'MovieData'));
ip.addRequired('MD_struct', @(x) isa(x,'MovieData'));
ip.addOptional('funParams',[],@isstruct);
ip.parse(MD_tracks, MD_struct, varargin{:});
%set default parameter if empty
parameter = ip.Results.funParams;
if isempty(parameter)
    parameter = PartitionAnalysisProcess.getDefaultParams(MD_tracks.outputDirectory_);
end
%% Initialization
%Easier variable access
channel_tracks = parameter.channel_tracks;
channel_struct = parameter.channel_struct;
outputDirectory = parameter.outputDirectory;
nControl = parameter.nControl;
%Get the indices of any previous partitioning processes from this function 
iPartProcess = MD_tracks.getProcessIndex('PartitionAnalysisProcess', 1, 0);
%If it doesn't exist already, create a new one
if isempty(iPartProcess)
    iPartProcess = numel(MD_tracks.processes_) + 1;
    MD_tracks.addProcess(PartitionAnalysisProcess(MD_tracks, MD_tracks.outputDirectory_));
end
partProcess = MD_tracks.processes_{iPartProcess};
%Check tracking process
iTrackProcess = MD_tracks.getProcessIndex('TrackingProcess');
if isempty(iTrackProcess)
    error('Tracking analysis has not been run! Please run tracking prior to partitioning analysis');
end
trackProcess = MD_tracks.processes_{iTrackProcess};
if ~trackProcess.checkChannelOutput(channel_tracks);
    error('Missing detection output! Please apply detection before running partitioning analysis');
end
%Check masking process
iMaskProcess = MD_struct.getProcessIndex('SubcellMaskProcess');
if isempty(iMaskProcess)
    error('Mask detected structure has not been run! Please run masking prior ro partitioning analysis');
end
maskProcess = MD_struct.processes_{iMaskProcess};
if ~maskProcess.checkChannelOutput(channel_struct);
    error('Missing masking output! Please apply detection before running partitioning analysis');
end
%Obtains data from tracking and masking processes
tracks = trackProcess.loadChannelOutput(1);
mask = maskProcess.loadChannelOutput(1);
%Checks if MD_tracks and MD_struct are compatible
%same imSize, and nFrames or nFrame of mask = 1
if ~all(MD_tracks.imSize_ == MD_struct.imSize_)
    error('The image sizes do not match');
end
if ~(MD_tracks.nFrames_ == MD_struct.nFrames_ || MD_struct.nFrames_ == 1)
    error('The movie lengths are not compatible');
end
%Set up the input directories
inFilePaths = cell(2,numel(MD_tracks.channels_));
inFilePaths{1, channel_tracks} = trackProcess.outFilePaths_{1, channel_tracks};
inFilePaths{2, channel_tracks} = maskProcess.outFilePaths_{1, channel_struct};
partProcess.setInFilePaths(inFilePaths);
%Set up the output file
outFilePaths = cell(1,numel(MD_tracks.channels_));
outFileName = ['Channel_' num2str(channel_tracks) '_partitioning_result'];
outFilePaths{1, channel_tracks} = [outputDirectory filesep outFileName '.mat'];
mkClrDir(outputDirectory);
partProcess.setOutFilePaths(outFilePaths);
%Easier variable access
%In image coordinates
xMax = MD_tracks.imSize_(2);
yMax = MD_tracks.imSize_(1);
isSingleFrame = MD_struct.nFrames_ == 1;
%% Partitioning Analysis
%get ROI mask
iProcMask = MD_tracks.getProcessIndex('ImportCellMaskProcess',1,0); %cell mask
if ~isempty(iProcMask)
    ROIMask = imread(fullfile(MD_tracks.processes_{iProcMask}.funParams_.OutputDirectory,'cellMask_channel_1.tif'));
else
    ROIMask = true(yMax, xMax);
end
%ProgressText
progressTextMultiple('Analyzing MD', nControl + 1);
%calls function that does partititoning analysis
partitionResult = trackPartitioning_StandAlone(tracks, mask, ROIMask, xMax, yMax, isSingleFrame, 'scrambleTracks', false); %#ok<NASGU>
progressTextMultiple();
%do randomized control
partitionControl = cell(1, nControl);
for iControl = 1:nControl
    partitionControl{iControl} = trackPartitioning_StandAlone(tracks, mask, ROIMask, xMax, yMax, isSingleFrame, 'scrambleTracks', true); %%#ok<NASGU>
    progressTextMultiple();
end
%partitionControl = vertcat(partitionControl_{:}); %#ok<NASGU>
%% Saving
save(outFilePaths{1,channel_tracks}, 'partitionResult', 'partitionControl');
MD_struct.save;
MD_tracks.save;
end

