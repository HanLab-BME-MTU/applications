function detectMoviekEBs(movieData,varargin)
%DETECTMOVIEKEBS detects EB signal at detected kinetochores
%
% Khuloud Jaqaman, Sebastien Besson, October 2012

%% Input

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn = ip.Results.paramsIn;

%Get the indices of any previous kinetochore-EB detection processes                                                                         
%If the process doesn't exist, create it
iProc = movieData.getProcessIndex('KEBDetectionProcess',1,0);
if isempty(iProc)
    iProc=numel(movieData.processes_)+1;
    movieData.addProcess(KEBDetectionProcess(movieData));
end
detProc = movieData.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(detProc,paramsIn);

%% --------------- Initialization ---------------%%

% Get grouping process index
iGroupProc = movieData.getProcessIndex('SisterGroupingProcess',1,1);

assert(~isempty(iGroupProc),['Sister grouping has not been run! '...
    'Please run sister grouping prior to k-EB detection!']);
groupProc=movieData.processes_{iGroupProc};

assert(all(groupProc.checkChannelOutput(setdiff([1 2],p.ChannelIndex))),...
    ['Missing sister grouping output ! Please apply sister grouping before ' ...
    'running  k-EB detection!']);

% Get tracking process index
iTrackProc = movieData.getProcessIndex('TrackingProcess',1,1);

assert(~isempty(iTrackProc),['Tracking has not been run! '...
    'Please run tracking prior to k-EB detection!'])
trackProc = movieData.processes_{iTrackProc};

assert(all(trackProc.checkChannelOutput(setdiff([1 2],p.ChannelIndex))),...
    ['Missing tracking output ! Please apply tracking before ' ...
    'running  k-EB detection!']);

% Set up the input directories (EB images + CENPA sisters + CENPA tracks)
inFilePaths = cell(3,numel(movieData.channels_));
for i = p.ChannelIndex
    inFilePaths{1,i} = movieData.getChannelPaths{i};
    inFilePaths{2,i} = groupProc.outFilePaths_{1,setdiff([1 2],i)};
    inFilePaths{3,i} = trackProc.outFilePaths_{1,setdiff([1 2],i)};
end
detProc.setInFilePaths(inFilePaths);

% Set up the output file
outFilePaths = cell(2,numel(movieData.channels_));
mkClrDir(p.OutputDirectory)
for i = p.ChannelIndex;
    outFilePaths{1,i} = [p.OutputDirectory filesep 'channel_' num2str(i) '.mat'];
    outFilePaths{2,i} = [p.OutputDirectory filesep 'filtered_images_for_channel_' num2str(i)];
    mkClrDir(outFilePaths{2,i})
end
detProc.setOutFilePaths(outFilePaths);

%% --------------- k-EB detection ---------------%%%

disp('Starting detecting k-EB signal ...');

%get image size and number of images
nImages = movieData.nFrames_;   
imSize = movieData.imSize_;

for i = p.ChannelIndex
    
    %read input
    tracks = trackProc.loadChannelOutput(setdiff([1 2],i));
    [sisterList,trackPairs] = groupProc.loadChannelOutput(setdiff([1 2],i));
    imageEB = zeros(imSize(1),imSize(2),nImages);
    for iImage = 1 : nImages
        imageEB(:,:,iImage) = movieData.channels_(i).loadImage(iImage);
    end
    
    %call function to detect kinetochore-EB signal
    sisterListEB = detectkEBs(sisterList,trackPairs,tracks,imageEB,...
        'radiusEB',p.radiusEB);
    stdList = NaN(nImages,1);
        
    %save output
    save(outFilePaths{1,i} ,'sisterListEB','stdList');
    
end

disp('Finished detecting k-EB signal!');

