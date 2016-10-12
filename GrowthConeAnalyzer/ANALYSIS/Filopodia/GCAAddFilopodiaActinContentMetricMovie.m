function [ output_args ] = GCAAddFilopodiaActinContentMetricMovie(movieData,varargin)
%GCAAddActinContentMetricMovie(movieData) 

% for now check movieData separately.
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end
%%Input check
ip = inputParser;

ip.CaseSensitive = false;

% PARAMETERS
defaultInDir = [movieData.outputDirectory_ filesep 'SegmentationPackage' ... 
    filesep 'StepsToReconstruct' filesep... 
    'VII_filopodiaBranch_fits' filesep 'Channel_1']; 
ip.addParameter('InputDirectory',defaultInDir,@(x) ischar(x)); 
ip.addParameter('ChannelIndex',1); 
ip.addParameter('ChannelIndexVeil',1); 

ip.addParameter('frames',[]); 

ip.parse(varargin{:});
p = ip.Results;
%% Initiate
load([ip.Results.InputDirectory filesep 'filoBranch.mat'])
load([movieData.outputDirectory_ filesep .... 
    'SegmentationPackage' filesep 'StepsToReconstruct' filesep ... 
    'III_veilStem_reconstruction' filesep 'Channel_' num2str(ip.Results.ChannelIndexVeil) filesep 'veilStem.mat']); 

% 
if isempty(ip.Results.frames)
    frames = 1:length(filoBranch)-1;
else
    frames = ip.Results.frames;
end

% find the background subtraction processes
idxProc = find(cellfun(@(x)  strcmpi(x.name_,'Background Subtraction'),movieData.processes_));
if length(idxProc)> 1;
    dispay('More than one background subtraction process found: Taking last');
    idxProc = idxProc(1);
   
end

if ~isempty(idxProc)
    imgDir = [ movieData.processes_{idxProc}.funParams_.OutputDirectory filesep 'background_subtracted_images_for_channel_' num2str(ip.Results.ChannelIndex)];
    imgFilenames = searchFiles('.tif',[],imgDir,0,'all',1);
     toAdd = 'BS';
else % use the channel data
    imgDir = [movieData.channels_(ip.Results.ChannelIndex).channelPath_]; 
    imgFilenames = searchFiles('.tif',[],imgDir,0,'all',1); 
    toAdd = []; 
end

for iFrame = 1:length(frames)
    
    % extract the veil
    veilMask = veilStem(frames(iFrame)).finalMask;
    % extract the img to feed into the function
    % changed 20160711 from img to background subtracted img 
    % img = double(imread([movieData.getChannelPaths{1} filesep movieData.getImageFileNames{1}{frames(iFrame)}]));
    img = double(imread(imgFilenames{frames(iFrame)}));  % load the background subtracted image 
  
    % extract the filo info to read into the function
    filoInfo = filoBranch(frames(iFrame)).filoInfo;
    % add the metric to the filo info - NOTE in the future might want to just
    % calculate automaticaly at the time of fitting to be more efficient.
    [filoInfo,normFactPerFrame] = GCAAddFilopodiaActinContentMetric(img,veilMask,filoInfo);
    filoBranch(frames(iFrame)).filoInfo = filoInfo;
    measC{frames(iFrame)} = normFactPerFrame;
end
% resave the values 
save([ip.Results.InputDirectory filesep 'filoBranch.mat'],'filoBranch','-v7.3'); 

%% save the normalization factor in the measurements folder 

% NOTE For final change to MEASUREMENT_EXTRACTION 
 expFolder = [movieData.outputDirectory_ filesep 'MEASUREMENT_EXTRACTION' filesep 'Descriptor' filesep 'GrowthCone' filesep 'ExpressionNormalization' toAdd]; 
% 
if ~isdir(expFolder) 
    mkdir(expFolder);
end 
 save([expFolder filesep 'meas_ExpressionNormalization.mat'],'measC'); 

end

