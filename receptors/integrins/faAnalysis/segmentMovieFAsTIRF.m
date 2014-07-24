function segmentMovieFAsTIRF(MD,thresholdMethod,methodValue,filterNoise,...
    filterBackground,minSize,plotRes)
%SEGMENTMOVIEFASTIRF ...
%
%SYNOPSIS segmentMovieFAsTIRF(MD,thresholdMethod,methodValue,filterNoise,...
%    filterBackground,minSize,plotRes)
%
%INPUT  MD    : The movieData object as output by the cell masking and
%               windowing software. Before calling this code,
%               thresholdMovie and refineMovieMask must have been
%               processed. Otherwise the code will crash.
%
%
%OUTPUT ...
%REMARKS The code will replace the refined masks with the real FA masks.
%
%Khuloud Jaqaman, November 2011

%% Input/Output

if nargin < 1
    error('segmentMovieFAsTIRF: Wrong number of input arguments');
end

%get image and analysis directories
imageDir = MD.channels_.channelPath_;
analysisDir = MD.movieDataPath_;

%get image and mask file listings
imageFileListing = dir([imageDir filesep '*.tif']);
if isempty(imageFileListing)
    imageFileListing = dir([imageDir filesep '*.tiff']);
end
refinedMasksDir = [analysisDir filesep 'refined_masks'];
refinedMasksDirFull = [refinedMasksDir filesep 'refined_masks_for_channel_1'];
refinedMaskFileListing = dir([refinedMasksDirFull filesep '*.tif']);
if isempty(refinedMaskFileListing)
    refinedMaskFileListing = dir([refinedMasksDirFull filesep '*.tiff']);
end

%get number of files
numFiles = length(imageFileListing);

%% FA segmentation

wtBar = waitbar(0,'Please wait, segmenting FAs ...'); 

for iFile = 1 : numFiles
    
    waitbar(iFile/numFiles,wtBar,'Please wait, segmenting FAs ...'); 
    
    %read image and mask
    image = double(imread(fullfile(imageDir,imageFileListing(iFile).name)));
    
    %call FA segmentation function
    maskFAs = segmentFAsTIRF(image,thresholdMethod,methodValue,filterNoise,...
        filterBackground,minSize,plotRes);
    
    %store new mask
    imwrite(maskFAs,fullfile(refinedMasksDirFull,refinedMaskFileListing(iFile).name),'tif');
    
end

close(wtBar)

%% ~~~ the end ~~~

