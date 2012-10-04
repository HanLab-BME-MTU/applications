function groupMovieSisters(movieData,varargin)
% Track features in a movie which has been processed by a detection method
%
% Sebastien Besson, 5/2011

%% Input
%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous tracking processes from this function                                                                              
iProc = movieData.getProcessIndex('SisterGroupingProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(SisterGroupingProcess(movieData));
end

groupProc = movieData.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(groupProc,paramsIn);

%% --------------- Initialization ---------------%%

% Check detection process first
iTrackProc =movieData.getProcessIndex('TrackingProcess',1,1);

assert(~isempty(iTrackProc),['Tracking has not been run! '...
    'Please run tracking prior to sister grouping!'])
trackProc = movieData.processes_{iTrackProc};

assert(all(trackProc.checkChannelOutput(p.ChannelIndex)),...
    ['Missing tracking output ! Please apply tracking before ' ...
    'running  sister grouping!']);
    
% Set up the input directories (input images)
inFilePaths = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex
    inFilePaths{1,i} = trackProc.outFilePaths_{1,i};
end
groupProc.setInFilePaths(inFilePaths);
    
% Set up the output file
outFilePaths = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex
    outFilePaths{1,i} = [p.OutputDirectory filesep 'channel_' num2str(i) '.mat'];
end
mkClrDir(p.OutputDirectory);
groupProc.setOutFilePaths(outFilePaths);

%% --------------- Displacement field calculation ---------------%%% 

disp('Starting grouping sisters...')

for i = p.ChannelIndex    
    tracks = trackProc.loadChannelOutput(i);
    
    [sisterList,trackPairs] = groupSisters(tracks,movieData.nFrames_,0,...
        'maxAngle', p.maxAngle, 'maxDist', p.maxDist,...
        'minOverlap', p.minOverlap, 'robust', p.robust);   %#ok<NASGU,ASGLU>
    
    % save each projData in its own directory
    save(outFilePaths{1,i},'sisterList','trackPairs')
end


disp('Finished grouping sisters!')
