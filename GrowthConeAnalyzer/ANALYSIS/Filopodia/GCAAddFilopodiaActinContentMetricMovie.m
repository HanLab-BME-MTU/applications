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
    filesep 'StepsToReconstruct'... 
    'VII_filopodiaBranch_fits']; 
ip.addParameter('InputDirectory',defaultInDir,@(x) ischar(x)); 

ip.parse(varargin{:});
p = ip.Results;
%% Initiate
load([ip.Results.InputDirectory filesep 'filoBranch.mat'])


for iFrame = 1:length(filoBranch)-1
load([movieData.outputDirectory_ filesep .... 
    'SegmentationPackage' filesep 'StepsToReconstruct' filesep ... 
    'III_veilStem_reconstruction' filesep 'Channel_1' filesep 'veilStem.mat']);    
    
% extract the veil
veilMask = veilStem(iFrame).finalMask;
% extract the img to feed into the function
img = double(imread([movieData.getChannelPaths{1} filesep movieData.getImageFileNames{1}{iFrame}])); 
% extract the filo info to read into the function 
filoInfo = filoBranch(iFrame).filoInfo;
% add the metric to the filo info - NOTE in the future might want to just
% calculate automaticaly at the time of fitting to be more efficient. 
[filoInfo,normFactPerFrame] = GCAAddFilopodiaActinContentMetric(img,veilMask,filoInfo); 
filoBranch(iFrame).filoInfo = filoInfo; 
paramC{iFrame} = normFactPerFrame; 
end
% resave the values 
save([ip.Results.InputDirectory filesep 'filoBranch.mat'],'filoBranch','-v7.3'); 

%% save the normalization factor in the measurements folder 

% NOTE For final change to MEASUREMENT_EXTRACTION 
% expFolder = ['PARAMETER_EXTRACTION' filespe 'Descriptor' filespe 'GrowthCone' filesep 'ExpressionNormalization']; 
% 
% if ~isdir(measurementFolder) 
%     mkdir(measurementFolder);
% end 
% save([expFolder filesep 'param_ExpressionNormalization.mat'],'paramC'); 






end

