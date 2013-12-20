function smDensity = particleDensityFromTracksInMask(tracksFinal,firstMaskFile,lengthMinMax)
%PARTICLEDENSITYFROMTRACKSINMASK returns single molecule density within cell mask
%
%SYNOPSIS smDensity = particleDensityFromTracksInMask(tracksFinal,firstMaskFile)
%
%INPUT  tracksFinal  : Output of trackCloseGapsKalman.
%       firstMaskFile: Name of first mask file, including full path.
%       lengthMinMax : Row vector with 2 entries indicating minimum and
%                      maximum length of a trajectory to include in
%                      analysis.
%                      Optional. Default: [5 99].
%
%OUTPUT smDensity: Average single molecule density per celle edge frame
%
%Khuloud Jaqaman, October 2013

%% Input

if nargin < 2
    disp('--particleDensityFromTracksInMask: Incorrect number of input arguments!');
    return
end

if nargin < 3 || isempty(lengthMinMax)
    lengthMinMax = [5 99];
end

%% Total cell area

%get names of mask files
outFileList = getFileStackNames(firstMaskFile);
numMaskFrames = length(outFileList);

%read all masks
mask1 = imread(outFileList{1});
maskStack = repmat(mask1,[1 1 numMaskFrames]);
for iMaskFrame = 2 : numMaskFrames
    maskStack(:,:,iMaskFrame) = imread(outFileList{iMaskFrame});        
end

%calculate total cell area
cellArea = sum(double(maskStack(:)));

%% Density

criteria.lifeTime.min = lengthMinMax(1);
criteria.lifeTime.max = lengthMinMax(2);
indxMinMax = chooseTracks(tracksFinal,criteria);

smDensity = length(indxMinMax) / cellArea;

%% ~~~ the end ~~~

