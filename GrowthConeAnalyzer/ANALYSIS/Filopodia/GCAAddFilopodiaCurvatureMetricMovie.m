function [ output_args ] = GCAAddFilopodiaCurvatureMetricMovie(movieData,varargin)
%GCAAddActinContentMetricMovie(movieData) 
% 

%% Input check
% for now check movieData separately.
if nargin < 1 || ~isa(movieData,'MovieData')
    error('The first input must be a valid MovieData object!')
end

ip = inputParser;

ip.CaseSensitive = false;

% PARAMETERS
defaultInDir = [movieData.outputDirectory_ filesep 'SegmentationPackage' ... 
    filesep 'StepsToReconstruct' filesep... 
    'VII_filopodiaBranch_fits' filesep 'Channel_1']; 
ip.addParameter('InputDirectory',defaultInDir,@(x) ischar(x)); 

ip.parse(varargin{:});

%% Initiate
load([ip.Results.InputDirectory filesep 'filoBranch.mat']); 

%% Wrap

for iFrame = 1:length(filoBranch)-1
% extract the veil
%veilMask = analInfo(iFrame).masks.neuriteEdge;
% extract the img to feed into the function
%img = double(imread([movieData.getChannelPaths{1} filesep movieData.getImageFileNames{1}{iFrame}])); 
% extract the filo info to read into the function 
filoInfo = filoBranch(iFrame).filoInfo;
% add the metric to the filo info - NOTE in the future might want to just
% calculate automaticaly at the time of fitting to be more efficient. 

%   img =   imread([MD.channelPath_ filesep MD.getImageFileNames{1}{iFrame}]); 
%   imshow(-img,[]); 
%   arrayfun(@(x) 

filoInfo = GCAAddFilopodiaCurvature(filoInfo); 
filoBranch(iFrame).filoInfo = filoInfo; 

end

%% Save all the information per a given movie in the filoBranch Structure  
save([ip.Results.InputDirectory filesep 'filoBranch.mat'],'filoBranch','-v7.3'); 





end

