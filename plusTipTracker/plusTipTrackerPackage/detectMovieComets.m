function detectMovieComets(movieData,varargin)
% detectMovieComets detect comets in a movie
%
% detectMovieComets 
%
% SYNOPSIS detectMovieComets(movieData,paramsIn)
%
% INPUT   
%   movieData - A MovieData object describing the movie to be processed
%
%   paramsIn - Structure with inputs for optional parameters. The
%   parameters should be stored as fields in the structure, with the field
%   names and possible values as described below
%
% OUTPUT   

% Sebastien Besson, Sep 2011 (last modified Sep 2011)

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous speckle detection processes                                                                     
iProc = movieData.getProcessIndex('CometDetectionProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(CometDetectionProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end
cometDetProc = movieData.processes_{iProc};
%Parse input, store in parameter structure
p = parseProcessParams(cometDetProc,paramsIn);

%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...','Name',cometDetProc.getName());
else 
    wtBar = -1;
end

% Reading various constants
imDirs  = movieData.getChannelPaths;
bitDepth = movieData.camBitdepth_;
nFrames = movieData.nFrames_;
maxIntensity =(2^bitDepth-1);
nChan =  numel(movieData.channels_); 

% Set up the input directories
inFilePaths = cell(1,nChan);
for j = p.ChannelIndex,  inFilePaths{1,j} = imDirs{j}; end
cometDetProc.setInFilePaths(inFilePaths);
    
% Set up the output directories
outFilePaths = cell(1,nChan);
mkClrDir(p.OutputDirectory);
for i = p.ChannelIndex;    
    %Create string for current directory
    outFilePaths{1,i} = [p.OutputDirectory filesep 'channel_' num2str(i) '.mat'];
end
cometDetProc.setOutFilePaths(outFilePaths);

%% --------------- Comet detection ---------------%%% 
disp('Starting detecting comets...')

logMsg = @(chan) ['Please wait, detecting comets for channel ' num2str(chan)];
timeMsg = @(t) ['\nEstimated time remaining: ' num2str(round(t)) 's'];
tic;

% Read region of interest
roiMask = movieData.getROIMask;
allMovieInfo(nFrames,nChan) = struct('xCoord', [], 'yCoord', [], 'amp', [],...
        'int',[],'ecc',[]);
    
for iChan= p.ChannelIndex    
    if ishandle(wtBar), waitbar(0,wtBar,'Loading image stack'); end    
    stack = movieData.channels_(iChan).loadImage(1:nFrames)/maxIntensity;
    
    
    % get difference of Gaussians image for each frame and standard deviation
    % of the cell background, stored in stdList
    stdList=nan(nFrames,1);
    filteredStack=zeros([movieData.imSize_ nFrames]);
    
    logMsg='Filtering images for comet detection';
    progressText(0,logMsg);
    if ishandle(wtBar), waitbar(0,wtBar,logMsg); end  
    for i = 1:nFrames
        if ~isempty(p.MaskProcessIndex)
            mask= roiMask(:,:,i) & movieData.processes_{p.MaskProcessIndex}.loadChannelOutput(iChan,i);
        else
            mask=roiMask(:,:,i);
        end
        % create kernels for gauss filtering
        blurKernelLow  = fspecial('gaussian', 21, p.sigma1);
        blurKernelHigh = fspecial('gaussian', 21, p.sigma2);
        
        Wlow = imfilter(double(mask), blurKernelLow);
        Whigh = imfilter(double(mask), blurKernelHigh);
        
        % use subfunction that calls imfilter to take care of edge effects
        % for now don't apply roiMask
        im=stack(:,:,i);
        im(~mask)=0;
        lowPass = imfilter(im, blurKernelLow)./Wlow;
        lowPass(~mask)=NaN;
        highPass = imfilter(im, blurKernelHigh)./Whigh;
        highPass(~mask)=NaN;
        filterDiff=lowPass-highPass;

        % Replace the imfilter
        % filterDiff=filterGauss2D(stack(:,:,i),p.sigma1) - filterGauss2D(stack(:,:,i),p.sigma2);
        
        % if there is a mask for each image file
        stdList(i)=nanstd(filterDiff(:));
        % (just removing not a numbers here from filterDiff so can take std)
        
        filteredStack(:,:,i)=filterDiff;
        progressText(i/nFrames,'Filtering images for comet detection');
        if ishandle(wtBar) && mod(i,5)==0, waitbar(i/nFrames,wtBar,logMsg); end  
    end
    
    meanStd = arrayfun(@(x) mean(stdList(max(1,x-1):min(nFrames,x+1))),1:nFrames);
    %         meanStd = smooth(stdList,3); % Discrepancy for endpoint
    
    % loop thru frames and detect
    logMsg='Detecting comets';
    progressText(0,logMsg);
    if ishandle(wtBar), waitbar(0,wtBar,logMsg); end  
    for i = 1:nFrames
        stepSize=p.multFactorStepSize*meanStd(i);
        thresh= p.multFactorThresh*meanStd(i);
        allMovieInfo(i,iChan) = detectComets(filteredStack(:,:,i),stepSize,thresh);
        progressText(i/nFrames,'Detecting comets');
        if ishandle(wtBar) && mod(i,5)==0, waitbar(i/nFrames,wtBar,logMsg); end
    end
    
    movieInfo=allMovieInfo(:,iChan); %#ok<NASGU>
    save(outFilePaths{1,iChan} ,'movieInfo','stdList','filteredStack');
end

if ishandle(wtBar), close(wtBar); end
disp('Finished detecting comets');
