function postProcessMovieComets(movieData,varargin)
% Track features in a movie which has been processed by a detection method
%
% Sebastien Besson, Feb 2012

%% Input
%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;


%Get the indices of any previous tracking processes from this function                                                                              
iProc = movieData.getProcessIndex('CometPostTrackingProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(CometPostTrackingProcess(movieData));                                                                                                 
end
postProc = movieData.processes_{iProc};

%Parse input, store in parameter structure
p = parseProcessParams(postProc,paramsIn);

%% --------------- Initialization ---------------%%

assert(~isempty(movieData.timeInterval_));
assert(~isempty(movieData.pixelSize_));

% Check detection process first
iTrackProc =movieData.getProcessIndex('TrackingProcess',1,1);

assert(~isempty(iTrackProc),['Tracking has not been run! '...
    'Please run tracking prior to post-processing!'])
trackProc = movieData.processes_{iTrackProc};

assert(all(trackProc.checkChannelOutput(p.ChannelIndex)),...
    ['Missing tracking output ! Please apply tracking before ' ...
    'running  post-processing!'])

iDetProc =movieData.getProcessIndex('CometDetectionProcess',1,1);
detProc=movieData.processes_{iDetProc};

assert(all(detProc.checkChannelOutput(p.ChannelIndex)),...
    ['Missing detection output ! Please apply detection before ' ...
        'running post-processing!']);
    
% Set up the input directories (input images)
inFilePaths = cell(2,numel(movieData.channels_));
for i = p.ChannelIndex
    inFilePaths{1,i} = trackProc.outFilePaths_{1,i};
    inFilePaths{2,i} = detProc.outFilePaths_{1,i};
end
postProc.setInFilePaths(inFilePaths);
    
% Set up the output file
outFilePaths = cell(2,numel(movieData.channels_));
for i = p.ChannelIndex
    outFilePaths{1,i} = [p.OutputDirectory filesep 'channel_' num2str(i) '.mat'];
    outFilePaths{2,i} = [p.OutputDirectory filesep 'channel_' num2str(i) '.txt'];
end
mkClrDir(p.OutputDirectory);
postProc.setOutFilePaths(outFilePaths);

%% --------------- Displacement field calculation ---------------%%% 

disp('Starting post-processing...')

for i = p.ChannelIndex
    movieInfo = detProc.loadChannelOutput(i);
    tracksFinal = trackProc.loadChannelOutput(i);
    
    projData.secPerFrame = movieData.timeInterval_;
    projData.pixSizeNm = movieData.pixelSize_;

    % figure out which frames were used in detection
    detExists=find(arrayfun(@(x) ~isempty(x.xCoord),movieInfo));
    sF=min(detExists); eF=max(detExists);
    
    % frame ranges for each step
    projData.detectionFrameRange=[sF eF];
    
    [projData,M]=postProcessMTTracks(projData,tracksFinal,movieInfo,[]);
    
    % save each projData in its own directory
    save(outFilePaths{1,i},'projData')
    
    % write out speed/lifetime/displacement distributions into a text file
    dlmwrite(outFilePaths{2,i}, M,...
        'precision', 3,'delimiter', '\t','newline', 'pc');
    
    % Write stats results into a text file
    statsFile = [p.OutputDirectory filesep 'Channel' num2str(i) '_Stats.txt'];
    statsData= struct2cell(projData.stats);
    statsName = fieldnames(projData.stats);
    fid=fopen(statsFile,'w+');
    for j=1:numel(statsName)
        fprintf(fid,'%s\t%g\n',statsName{j},statsData{j});
    end
    fclose(fid);
    
    if p.makeHist==1
        plusTipMakeHistograms(M,[p.OutputDirectory filesep 'histograms'])
%         plusTipPlotTrackAngles(runInfo,[p.OutputDirectory filesep 'histograms']);
    end
end


disp('Finished post-processing comets!')
