function detectMovieThreshLocMax(movieData, varargin)
%DETECTMOVIETHRESHLOCMAX compiles detection data from movie frames

%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.parse(movieData,varargin{:});
paramsIn=ip.Results.paramsIn;

%Get the indices of any previous speckle detection processes                                                                     
iProc = movieData.getProcessIndex('ThreshLocMaxProcess',1,0);

%If the process doesn't exist, create it
if isempty(iProc)
    iProc = numel(movieData.processes_)+1;
    movieData.addProcess(ThreshLocMaxProcess(movieData,...
        movieData.outputDirectory_));                                                                                                 
end
threshLocMaxProc = movieData.processes_{iProc};
%Parse input, store in parameter structure
p = parseProcessParams(threshLocMaxProc,paramsIn);

%% Detection Parameters
 thresholdMethod=p.detectionParam.thresholdMethod;
 methodValue= p.detectionParam.methodValue;
 filterNoise=p.detectionParam.filterNoise; %default 1
 filterBackground=p.detectionParam.filterBackground; %default 10
 minSize=p.detectionParam.minSize;
 maxSize = p.detectionParam.maxSize;
 alphaLocMax=p.detectionParam.alphaLocMax;
 alphaSubRes=p.detectionParam.alphaSubRes;
 psfSigma = p.detectionParam.psfSigma;
 
 % Was masking done?
%  iProc = movieData.getProcessIndex('SegmentationProcess',1,0);
%  if ~isempty(iProc)
%      pMask = parseProcessParams(movieData.processes_{iProc},[]);
%      maskDir = movieData.processes_{iProc}.outFilePaths_(pMask.ChannelIndex);
%      aux = dir([maskDir{1} ,filesep,'*.tif']);
%      mask=1;
%  else
    mask=[];
%  end
%% --------------- Initialization ---------------%%

nChan=numel(movieData.channels_);
% Set up the input directories
inFilePaths = cell(1,nChan);
for i = p.ChannelIndex
    inFilePaths{1,i} = movieData.getChannelPaths{i};
end
threshLocMaxProc.setInFilePaths(inFilePaths);
    
% Set up the output directories
outFilePaths = cell(1,nChan);
saveResults(nChan,1)=struct();
dName = 'detections_for_channel_';
for i = p.ChannelIndex;    
    currDir = [p.OutputDirectory filesep dName num2str(i)];
    saveResults(i).dir = currDir ;
    saveResults(i).filename = ['Channel_' num2str(i) '_detection_result.mat'];
    %Create string for current directory
    outFilePaths{1,i} = [saveResults(i).dir filesep saveResults(i).filename ];
    threshLocMaxProc.setOutFilePaths(outFilePaths{1,i},i);
    mkClrDir(currDir);
end
list = p.lastImageNum;
% list = p.lastImageNum-p.firstImageNum;
movieInfo = repmat(struct('xCoord',[],'yCoord',[],'zCoord',[],'amp',[],'size',[],'ecc',[]),list,1); 
pixelRef = cell(list,1);
pixelIndx = zeros(movieData.imSize_(1),movieData.imSize_(2),list);
%% --------------- Sub-resolution object detection ---------------%%% 
disp('Starting detecting diffraction-limited objects...')

for i = p.ChannelIndex
    disp(['Please wait, detecting objects for channel ' num2str(i)])
    disp(inFilePaths{1,i});
    disp('Results will be saved under:')
    disp(outFilePaths{1,i});
    
    % Retrieve information about the images
    if movieData.isOmero() || movieData.isBF()
        movieParam.channel = movieData.channels_(i);
    else
        [~, base, digits4Enum, ~] = getFilenameBody(movieData.getImageFileNames{i}{1});
        digits4Enum = length(digits4Enum);
        
        movieParam.imageDir = [inFilePaths{1,i} filesep];
        movieParam.filenameBase = base;
        movieParam.digits4Enum = digits4Enum;
    end    
    movieParam.firstImageNum=p.firstImageNum;
    movieParam.lastImageNum=p.lastImageNum;
    j=1;
    maskThresh = [];
    for k = movieParam.firstImageNum:movieParam.lastImageNum
        %Load Image
        I = movieData.channels_(i).loadImage(k); %This seems dumb, changing(p.ChannelIndex(i)) to i
        if isempty(mask)
        else
            mask = imread([maskDir{1},'/',aux(k).name]);
        end
        %Run Detection    
        [detectedFeatures,pixelPos,maskBlobs] = detectThreshLocMax(I,...
        thresholdMethod,methodValue,filterNoise,filterBackground,minSize,...
        alphaLocMax,[],mask,maxSize,maskThresh,alphaSubRes,psfSigma);
        maskThresh = maskBlobs; %If masking fails
        movieInfo(k) = detectedFeatures;
        pixelRef{k} = pixelPos;
        test = cell2mat(pixelPos);
        pixelMap = zeros(movieData.imSize_);
        pixelMap(test) = 1;
        pixelIndx(:,:,k) = pixelMap;
        j = j+1; 
    end
    save(strcat(saveResults(i).dir,'/',saveResults(i).filename),'movieInfo','pixelIndx','pixelRef');
end

disp('Finished detecting diffraction-limited objects...')
end