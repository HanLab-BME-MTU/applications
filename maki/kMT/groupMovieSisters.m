function groupMovieSisters(movieData,varargin)
%GROUPMOVIESISTERS groups sisters in a movie which has been processed by a tracking method
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

%Get the indices of any previous sister grouping processes                                                                         
%If the process doesn't exist, create it
iProc = movieData.getProcessIndex('SisterGroupingProcess',1,0);
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(SisterGroupingProcess(movieData));
end
groupProc = movieData.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(groupProc,paramsIn);

%% --------------- Initialization ---------------%%

% Check tracking process
iTrackProc = movieData.getProcessIndex('TrackingProcess',1,1);

assert(~isempty(iTrackProc),['Tracking has not been run! '...
    'Please run tracking prior to sister grouping!'])
trackProc = movieData.processes_{iTrackProc};

assert(all(trackProc.checkChannelOutput(p.ChannelIndex)),...
    ['Missing tracking output ! Please apply tracking before ' ...
    'running  sister grouping!']);

% Check spindle pole process if necessary
if p.useAlignment || p.associateSis2Pole
    
    iPoleProc = movieData.getProcessIndex('SpindlePolesEBProcess',1,1);
    
    assert(~isempty(iPoleProc),['Spindle poles have not been detected! '...
        'Please run spindle pole detection prior to sister grouping!'])
    poleProc = movieData.processes_{iPoleProc};
    
    assert(all(poleProc.checkChannelOutput(setdiff([1 2],p.ChannelIndex))),...
        ['Missing spindle pole output! Please detect spindle poles before ' ...
        'running sister grouping!']);
    
end

% Set up the input directories (tracks+spindle poles)
inFilePaths = cell(1+(p.useAlignment||p.associateSis2Pole),numel(movieData.channels_));
for i = p.ChannelIndex
    inFilePaths{1,i} = trackProc.outFilePaths_{1,i};
    if p.useAlignment || p.associateSis2Pole
        inFilePaths{2,i} = poleProc.outFilePaths_{1,setdiff([1 2],i)};
    end
end
groupProc.setInFilePaths(inFilePaths);
    
% Set up the output file
outFilePaths = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex
    outFilePaths{1,i} = [p.OutputDirectory filesep 'channel_' num2str(i) '.mat'];
end
mkClrDir(p.OutputDirectory);
groupProc.setOutFilePaths(outFilePaths);

%% --------------- Sister grouping ---------------%%% 

disp('Grouping sisters...')

for i = p.ChannelIndex    
    
    tracks = trackProc.loadChannelOutput(i);
    if p.useAlignment || p.associateSis2Pole
        poleInfo = poleProc.loadChannelOutput(setdiff([1 2],i));
    else
        poleInfo = [];
    end
    
    [sisterList,trackPairs] = groupSisters(tracks,movieData.nFrames_,poleInfo,0,...
        'maxAngle', p.maxAngle, 'maxDist', p.maxDist,...
        'minOverlap', p.minOverlap, 'useAlignment', p.useAlignment, ...
        'robust', p.robust,'associateSis2Pole',p.associateSis2Pole);   %#ok<NASGU,ASGLU>
    
    % save each projData in its own directory
    save(outFilePaths{1,i},'sisterList','trackPairs')
    
end

disp('Done')
