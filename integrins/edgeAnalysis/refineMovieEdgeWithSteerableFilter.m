function refineMovieEdgeWithSteerableFilter(MD,doPlot)
% function refineMovieEdgeWithSteerableFilter(imageDir,analysisDir,doPlot)
%REFINEMOVIEEDGEWITHSTEERABLEFILTER replaces simple masks with masks refined using steerable line filtering
%
%SYNOPSIS refineMovieEdgeWithSteerableFilter(imageDir,analysisDir)
%
%INPUT  MD    : The movieData object as output by the cell masking and
%               windowing software. Before calling this code,
%               thresholdMovie and refineMovieMask must have been
%               processed. Otherwise the code will crash.
%       doPlot: 1 to plot refined masks, 0 otherwise. Refined masks shown
%               in green, original masks shown in blue. Note that this
%               refinement comes on top of the "refineMovieMask"
%               refinement.
%
%OUTPUT No output variables. HOWEVER: the code will copy the original
%       masks and refined masks into new directories called masks_OLD and
%       refined_masks_OLD and will replace the old masks with the ones it
%       produces. After this, one should run refineMovieMasks one more time
%       to get back on track and use the rest of the windowing functions.
%
%Khuloud Jaqaman, November 2011

%% Input/Output

if nargin ~= 2
    error('refineMovieEdgeWithSteerableFilter: Wrong number of input arguments');
end

%get image and analysis directories
imageDir = MD.channels_.channelPath_;
analysisDir = MD.movieDataPath_;

%make new directories and copy old masks to them
masksDir = [analysisDir filesep 'masks'];
masksOldDir = [analysisDir filesep 'masks_OLD'];
mkdir(masksOldDir);
copyfile([masksDir filesep '*'],masksOldDir);
refinedMasksDir = [analysisDir filesep 'refined_masks'];
refinedMasksOldDir = [analysisDir filesep 'refined_masks_OLD'];
mkdir(refinedMasksOldDir);
copyfile([refinedMasksDir filesep '*'],refinedMasksOldDir);

%get image and mask file listings
imageFileListing = dir(imageDir);
masksDirFull = [masksDir filesep 'masks_for_channel_1'];
maskFileListing = dir(masksDirFull);
refinedMasksDirFull = [refinedMasksDir filesep 'refined_masks_for_channel_1'];
refinedMaskFileListing = dir(refinedMasksDirFull);

%remove all listings that are not tif files
imageFileListing = keepOnlyTiffFiles(imageFileListing);
maskFileListing = keepOnlyTiffFiles(maskFileListing);
refinedMaskFileListing = keepOnlyTiffFiles(refinedMaskFileListing);

%get number of files
numFiles = length(imageFileListing);

%% Mask refinement

%refine masks using steerable line filter
for iFile = 1 : numFiles
    
    %read image and mask
    image = double(imread(fullfile(imageDir,imageFileListing(iFile).name)));
    mask0 = double(imread(fullfile(refinedMasksDirFull,refinedMaskFileListing(iFile).name)));
    
    %call refinement function
    mask = refineEdgeWithSteerableFilter(mask0,image,doPlot);
    
    %store new mask
    imwrite(mask,fullfile(masksDirFull,maskFileListing(iFile).name),'tif');
    
end

%% Sub-function

function fileListing = keepOnlyTiffFiles(fileListing)

numFiles = length(fileListing);
goodFile = zeros(numFiles,1);
for iFile = 1 : numFiles
    fileName = fileListing(iFile).name;
    if (length(fileName) > 4) && ((strcmp(fileName(end-2:end),'tif')) || ...
            (strcmp(fileName(end-3:end),'tiff')))
        goodFile(iFile) = 1;
    end
end
fileListing = fileListing(goodFile==1);




%% ~~~ the end ~~~

