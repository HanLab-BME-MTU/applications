function detectMovieSpindlePolesEB(movieData,varargin)
%DETECTMOVIESPINDLEPOLESEB detects the spindle poles from EB images by calling detectSpindlePolesEB
%
%Khuloud Jaqaman, 10/2012

%% Input
%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous spindle pole processes from this function                                                                              
iProc = movieData.getProcessIndex('SpindlePolesEBProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(SpindlePolesEBProcess(movieData));
end

axisProc = movieData.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(axisProc,paramsIn);

%% --------------- Initialization ---------------%%

% Set up the input directories (EB images)
inFilePaths = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex
    inFilePaths{1,i} = movieData.getChannelPaths{i};
end
axisProc.setInFilePaths(inFilePaths);

% Set up the output file
outFilePaths = cell(1,numel(movieData.channels_));
for i = p.ChannelIndex
    outFilePaths{1,i} = [p.OutputDirectory filesep 'channel_' num2str(i) '.mat'];
end
mkClrDir(p.OutputDirectory);
axisProc.setOutFilePaths(outFilePaths);

%% --------------- Spindle pole detection ---------------%%% 

disp('Detecting spindle poles...')

%get image size and number of images
nImages = movieData.nFrames_;   
imSize = movieData.imSize_;

for i = p.ChannelIndex
    
    %read images
    imageEB = zeros(imSize(1),imSize(2),nImages);
    for iImage = 1 : nImages
        imageEB(:,:,iImage) = movieData.channels_(i).loadImage(iImage);
    end
    
    %call function to detect spindle poles
    poleInfo = detectSpindlePolesEB(imageEB,'doPlot',p.doPlot,'numPoles',p.numPoles); %#ok<NASGU>
    
    %save each projData in its own directory
    save(outFilePaths{1,i},'poleInfo')
    
end

disp('Done')
