function groupMovieTracks(movieData,varargin)
% groupMovieTracks groups tracks using maximum 
%
% SYNOPSIS groupMovieTracks(movieData,paramsIn)
%
% INPUT   
%   movieData - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%
% OUTPUT   

% Sebastien Besson, Jun 2012

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous speckle detection processes                                                                     
iProc = movieData.getProcessIndex('TrackGroupingProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(TrackGroupingProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end
groupProc = movieData.processes_{iProc};
%Parse input, store in parameter structure
p = parseProcessParams(groupProc,paramsIn);

%% --------------- Initialization ---------------%%


%Get the detection process 
if ~isempty(p.DetProcessIndex)
    iDetProc = movieData.getProcessIndex('AnisoGaussianDetectionProcess',1,0);
    assert(~isempty(iDetProc),'Detection must be run prior to track grouping');
    p.DetProcessIndex = iDetProc;
end
detProc = movieData.processes_{p.DetProcessIndex};

%Get the tracking process                                                                    
if ~isempty(p.TrackProcessIndex)
    iTrackProc = movieData.getProcessIndex('TrackingProcess',1,0);
    assert(~isempty(iTrackProc),'Tracking must be run prior to track grouping');
    p.TrackProcessIndex = iTrackProc;
end
trackProc = movieData.processes_{p.TrackProcessIndex};

nChan=numel(movieData.channels_);
% Set up the input directories
inFilePaths = cell(2,nChan);
for i = p.ChannelIndex
    inFilePaths{1,i} = detProc.outFilePaths_{1,i};
    inFilePaths{2,i} = trackProc.outFilePaths_{1,i};
end
groupProc.setInFilePaths(inFilePaths);
    
% Set up the output directories
outFilePaths = cell(2,nChan);
mkClrDir(p.OutputDirectory);
for i = p.ChannelIndex;    
    %Create string for current directory
    outputPath = fullfile(p.OutputDirectory,['Channel_' num2str(i)]);
    if ~isdir(outputPath), mkdir(outputPath); end
    outFilePaths{1,i} = fullfile(outputPath, 'ClassifiedTracks.mat');
    outFilePaths{2,i} = fullfile(outputPath, 'ClassifiedSegments.mat');
    outFilePaths{3,i} = fullfile(outputPath, 'transformed_masks');
    mkClrDir(outFilePaths{3,i});
end
groupProc.setOutFilePaths(outFilePaths);

%% --------------- Sub-resolution object detection ---------------%%% 
disp('Starting grouping tracks...')

detParams = detProc.funParams_;
roiMask = movieData.getROIMask();
nFrames = movieData.nFrames_;
for i = p.ChannelIndex
    % Load input
    featuresInfo = detProc.loadChannelOutput(i);
    tracksFinal = trackProc.loadChannelOutput(i);

    fString = strcat('%0',num2str(ceil(log10(nFrames)+1)),'.f');

    for j=1:nFrames
        mask = movieData.processes_{p.MaskProcessIndex}.loadChannelOutput(p.ChannelIndex,j);
        mask = logical(mask) & roiMask(:,:,j);
        distToEdge = double(bwdist(~mask)); %#ok<NASGU>
    
        save(fullfile(outFilePaths{3,i},['distanceTransform_' num2str(j,fString) '.mat']),...
            'distToEdge')
    end

    
    getMoviePairTracks1(featuresInfo, tracksFinal, detParams.psfSigma, detParams.kSigma,...
        movieData.nFrames_, movieData.imSize_, movieData.pixelSize_, ...
        outFilePaths{3,i}, fullfile(p.OutputDirectory,['Channel_' num2str(i)]), p.minLifetime,...
        p.maxDistance, p.minOverlap, p.bandWidth, p. minDistance, ...
        p.alpha)
end

disp('Finished detecting objects...')

