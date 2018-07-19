function makiCheapVis(dataStruct,isCroppedMovie)
%MAKICHEAPVIS is a cheap way of visualizing spots for the maki project
%
% makiCheapVis creates a movie, placing Gaussians at the detected
%              coordinates, and writes it out as a TIFF-series into a
%              directory called "MOVIENAME_spotMovie". The directory is a
%              subdirectory to the dataFilePath. These series can be
%              loaded as an additional channel into your favorite 3D
%              visualization software to show the detection result.
%
% SYNOPSIS: makiCheapVis(dataStruct)
%
% INPUT dataStruct: maki data structure
%
% OUTPUT
%
% REMARKS Movies can be very big. Thereore, makiCheapVis generates and
%         each frame individually.
%
% created with MATLAB ver.: 7.5.0.342 (R2007b) on Windows_NT
%
% created by: Jonas Dorn
% DATE: 07-Mar-2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2 || isempty(isCroppedMovie)
    isCroppedMovie = false;
end

% get info from dataStruct. Read from movie header because of cropping.
if isobject(dataStruct)
    if ~isCroppedMovie
        imSize = dataStruct.imageData.imageProperties.imageSize(1:3);
    else
        imSize = dataStruct.imageSize(1:3);
    end
    nTimepoints = dataStruct.nTimepoints;
    imageName = dataStruct.identifier;
    sigma = dataStruct.dataProperties.FT_SIGMA;
    crop = dataStruct.dataProperties.imageData.cropInfo.crop(:,1:3);
else
    imSize = [dataStruct.movieHeader.numRows,dataStruct.movieHeader.numCols,dataStruct.movieHeader.numZSlices];
    nTimepoints = dataStruct.movieHeader.numTimepoints;
    sigma = dataStruct.dataProperties.FT_SIGMA;
    imageName = dataStruct.projectName;

    if isempty(dataStruct.dataProperties.crop)
        dataStruct.dataProperties.crop = zeros(2,3);
    end
    crop = dataStruct.dataProperties.crop(:,1:3);
end
% if the movie is cropped, we will not shift the coordinates
if isCroppedMovie
    isCrop = false(1,3);
else
    isCrop = any(crop,1);
end
crop(1,~isCrop) = 1;
crop(2,~isCrop) = imSize(~isCrop);



% create subdirectory
saveDirectory = fullfile(dataStruct.dataFilePath, ...
    sprintf('%s_spotMovie',imageName));

% -- here, I should check whether the directory exists already
mkdir(saveDirectory);

% set dynamic range. The images will be stored as 16 bit tiffs. Before
% making the movie, it cannot be known what the maximum intensity turns out
% to be. However, the code does not support superresolution yet, therefore
% the overlaps between signals are going to be minimal.
% Make the maximum spot intensity = 3500 greyvalues, and set the zero at 0.
allAmp = catStruct(1,'dataStruct.initCoord.amp(:,1)');
maxInt = max(allAmp);
minInt = 0;
nDigitsT=floor(log10(nTimepoints))+1;
nDigitsZ=floor(log10(imSize(3)))+1;
spstr=sprintf('%%s_z%%.%dd_t%%.%dd.tif',nDigitsZ,nDigitsT);

% loop and create the GaussImages
progressText(0,'Creating tiff series')
for t = 1:nTimepoints
    % create image
    if ~isempty(dataStruct.initCoord(t))
        coords = dataStruct.initCoord(t).allCoordPix(:,1:3);
        coords = coords + repmat(crop(1,:)-1,size(coords,1),1);
        image = putSpotsND(coords,...
            dataStruct.initCoord(t).initAmp(:,1),...
            imSize, sigma);
        % convert to uint16
        image = uint16((image-minInt)/maxInt*3500);
        % write image to file
        for z = 1:imSize(3)
            imwrite(image(:,:,z),fullfile(saveDirectory,sprintf(spstr,imageName,z,t)));
        end
        progressText(t/nTimepoints)
    end
end