function [imageData, dataProperties] = read3dFISH(varargin)
%READ3DFISH Reads microscope images and give imagedata 3D stach with data
%properties
%   Input(optional): pathName

%   Output: imageData: Multichannel 3D stack (for findNuclei.m)
%           dataProperties (for singleNucleusSpotDetection.m)
% 
% 03/2016 Ning
% 
% p = inputparser;
% p.addOptional('history', '', @ischar);

autoloadBioFormats = 1;
% load the Bio-Formats library into the MATLAB environment
status = bfCheckJavaPath(autoloadBioFormats);
assert(status, ['Missing Bio-Formats library. Either add bioformats_package.jar '...
    'to the static Java path or add it to the Matlab path.']);

if nargin == 0
    [fileName, pathName] = uigetfile(bfGetFileExtensions, 'Select the file contaning multichannel 3D stack');
else
    pathName = varargin{1};
    [fileName, pathName] = uigetfile(bfGetFileExtensions, 'Select the file contaning multichannel 3D stack', pathName);
end

% Try to store and retrieve path history in a better way
save pathName pathName

dataFilePath = fullfile( pathName, fileName );
MD = MovieData.load(dataFilePath);
% reader = MD.getReader();

% Define dataProperties parameters
% Pixelsize is in um
dataProperties.PIXELSIZE_XY = MD.pixelSize_/1000;
dataProperties.PIXELSIZE_Z = MD.pixelSizeZ_/1000;
dataProperties.imSize = MD.imSize_;
dataProperties.nDepth = MD.zSize_;
dataProperties.NA = MD.numAperture_;

lenseType = input('Enter lense type (air, water or oil) > ', 's');
switch lenseType
    case 'air'
        dataProperties.refractiveIndex = 1;
    case 'water'
        dataProperties.refractiveIndex = 1.33;
    case 'oil'
        dataProperties.refractiveIndex = 1.51;
end
% sigmaCorrection defined by default
dataProperties.sigmaCorrection=[1 1];



dapiChaIndex = input('Enter the index number (1, 2, 3...) for dapi channel > ');
dapi = MD.getChannel(dapiChaIndex);
% Get frame size of a single channel then load 3D stack for all frames
nFrameDapi = dapi.getReader().getSizeT;
imageData.dapi = dapi.loadStack(nFrameDapi);

% Figure out which wavelength (excitation or emission) is used for psf size
% calculation.
% PSF size is in um
dataProperties.dapiWVL = dapi.emissionWavelength_/1000;


greenChaIndex = input('Enter the index number (1, 2, 3...) for green channel > ');
green = MD.getChannel(greenChaIndex);
% Get frame size of a single channel then load 3D stack for all frames
nFrameGreen = green.getReader().getSizeT;
imageData.green = green.loadStack(nFrameGreen);

dataProperties.greenWVL = green.emissionWavelength_/1000;


redChaIndex = input('Enter the index number (1, 2, 3...) for red channel > ');
red = MD.getChannel(redChaIndex);
% Get frame size of a single channel then load 3D stack for all frames
nFrameRed = red.getReader().getSizeT;
imageData.red = red.loadStack(nFrameRed);

dataProperties.redWVL = red.emissionWavelength_/1000;

clear MD dapi green red

end

