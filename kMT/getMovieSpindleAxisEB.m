function getMovieSpindleAxisEB(movieData,varargin)
%GETMOVIESPINDLEAXISEB determines spindle axis from EB images
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

%Get the indices of any previous spindle axis processes from this function                                                                              
iProc = movieData.getProcessIndex('SpindleAxisEBProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(SpindleAxisEBProcess(movieData));
end

axisProc = movieData.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(axisProc,paramsIn);

%% --------------- Initialization ---------------%%

%MAYBE NEEDED IN THE FUTURE
%
% % Check detection process first
% iDetProc = movieData.getProcessIndex('DetectionProcess',1,1);
% 
% assert(~isempty(iDetProc),['Detection has not been run! '...
%     'Please run detection prior to spindle axis estimation!'])
% detProc = movieData.processes_{iDetProc};
% 
% assert(all(detProc.checkChannelOutput(p.ChannelIndex)),...
%     ['Missing detection output ! Please apply detection before ' ...
%     'running spindle axis estimation!']);
    
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

%% --------------- Spindle axis estimation ---------------%%% 

disp('Estimating spindle axis...')

%get image size and number of images
nImages = movieData.nFrames_;   
imSize = movieData.imSize_;

for i = p.ChannelIndex
    
    %read images
    imageEB = zeros(imSize(1),imSize(2),nImages);
    for iImage = 1 : nImages
        imageEB(:,:,iImage) = movieData.channels_(i).loadImage(iImage);
    end
    
    %call function to estimate spindle axis
    [spindleAxisVec,poleInfo] = getSpindleAxisEB(imageEB,'doPlot',p.doPlot); %#ok<ASGLU,NASGU>
    
    %save each projData in its own directory
    save(outFilePaths{1,i},'spindleAxisVec','poleInfo')
    
end

disp('Done')
