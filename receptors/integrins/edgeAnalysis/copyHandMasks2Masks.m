function copyHandMasks2Masks(MD)
%copyHandMasks2Masks 
%
%SYNOPSIS copyHandMasks2Masks(MD)
%
%INPUT  MD    : The movieData object as output by the cell masking and
%               windowing software. Before calling this code,
%               thresholdMovie and refineMovieMask must have been
%               processed. Otherwise the code will crash.
%
%OUTPUT 
%
%REMARKS The code will copy the original masks and refined masks into new
%        directories called masks_OLD2 and refined_masks_OLD2 and will
%        replace the old masks with those in groundTruthChannel1. After
%        this, one should run refineMovieMasks one more time to get back on
%        track and use the rest of the windowing functions.
%
%Khuloud Jaqaman, October 2012

%% Input/Output

if nargin < 1
    error('copyHandMasks2Masks: Wrong number of input arguments');
end

%get analysis directory
analysisDir = [MD.movieDataPath_ filesep 'SegmentationPackage'];

%make new directories and copy old masks to them
masksDir = [analysisDir filesep 'masks'];
masksOldDir = [analysisDir filesep 'masks_OLD2'];
mkdir(masksOldDir);
copyfile([masksDir filesep '*'],masksOldDir,'f');
refinedMasksDir = [analysisDir filesep 'refined_masks'];
refinedMasksOldDir = [analysisDir filesep 'refined_masks_OLD2'];
mkdir(refinedMasksOldDir);
copyfile([refinedMasksDir filesep '*'],refinedMasksOldDir,'f');

%get mask file listings
masksDirFull = [masksDir filesep 'masks_for_channel_1'];
maskFileListing = dir([masksDirFull filesep '*.tif']);
if isempty(maskFileListing)
    maskFileListing = dir([masksDirFull filesep '*.tiff']);
end
if isempty(maskFileListing)
    maskFileListing = dir([masksDirFull filesep '*.TIF']);
end

%get manually-curated masks
manualMasksDirFull = [analysisDir filesep 'groundTruthChannel1'];
manualMaskFileListing = imDir(manualMasksDirFull,0);

%get number of files
numFiles = length(maskFileListing);

%copy manual masks to refined masks directory
for iFile = 1 : numFiles
    copyfile(fullfile(manualMasksDirFull,manualMaskFileListing(iFile).name),fullfile(masksDirFull,maskFileListing(iFile).name),'f');
end

%% ~~~ the end ~~~

