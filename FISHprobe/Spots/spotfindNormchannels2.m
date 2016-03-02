function spots = spotfindNormchannels2(movieName,Cha)
startTime = clock;
%SPOTFINDTEST locates fluorescent tags in 3D data with modified
%dataProperties
%
% SYNOPSIS cordVar = spotfind(img)
%
% INPUT movieName   :   input movie data name
%       channelNum  :   movie channel number
%       MAXSPOTS    :   dataProperties.MAXSPOTS (suggest 5)
%       filterSize  :   index factor for FILTERPRM (suggest 4)
%
% OUTPUT spots : nTimepoints-by-1 structure with fields
%                .sp  structure with fields
%                   .cord coordinates
%                   .mnint spottiness
%                .COM center of mass of image

% 04/15 Ning
% Take movie data and process to get 3D matrix
Path = strcat('', movieName);
MD = MovieData.load(Path);
xSize = MD.imSize_(1);
ySize = MD.imSize_(2);
nDepth = MD.zSize_;
frameN = 1;

zStack3D = zeros(xSize, ySize, nDepth);
for i = 1:nDepth
    zStack3D(:,:,i) = MD.getChannel(Cha).loadImage(frameN,i);
end
zStack3D=100*zStack3D/max(zStack3D(:));

% ============================================================
% Define dataProperties parameters
% Generate dataProperties.FILTERPRM
if isempty(MD.numAperture_)
    dataProperties.NA = input('Enter Numerical Aperture > ');
else
    dataProperties.NA=MD.numAperture_;
end

if isempty(MD.channels_(1,Cha).emissionWavelength_)
    dataProperties.WVL = input('Enter Emission Wavelenth in um > ');
else
    dataProperties.WVL = MD.channels_(1,Cha).emissionWavelength_/1000;
end

dataProperties.sigmaCorrection=[1 1];

% Pixel size in um
% XY pixel size 10-fold calibration for nikon widefield
% dataProperties.PIXELSIZE_XY = MD.pixelSize_/1000*10;
dataProperties.PIXELSIZE_XY = MD.pixelSize_/1000;
dataProperties.PIXELSIZE_Z = MD.pixelSizeZ_/1000;

% Gaussian filter sigma 1-3 and size 4-6 for x y z respectively
% odd mask size with +/- 4*sigma is used
% Gussian filter with single nucleus size
nulDiameter = 15;
% nulDiameter = imput('Enter single nucleus diameter in um > ');
% Pixel size covered by a single nucleus
nucleus_xy = nulDiameter/dataProperties.PIXELSIZE_XY;
nucleus_z = nulDiameter/dataProperties.PIXELSIZE_Z;
% nucleus_xy = 87;
% nucleus_z = nucleus_xy/2200*nDepth;
fSze_xy = roundOddOrEven(nucleus_xy,'odd','inf');
fSze_z = roundOddOrEven(nucleus_z,'odd','inf');
gausspar = [fSze_xy/4,fSze_xy/4,fSze_z/4,fSze_xy,fSze_xy,fSze_z];
fImg = filtermovie(zStack3D,gausspar,0);


% Copied from defaultDataProperties.m
% calcFilterParms generates psf width, which is used as sigma for gaussian
% filter
% Refractive index: air 1; water 1.33 oil 1.51

% Code borrowed from 
% /home2/nzhang/matlab/applications/FISHprobe/Spots/detect/defaultDataProperties.m

% /home2/nzhang/matlab/common/mathfun/psfModels/calcFilterParms.m cauculate
% psf size as patchsize
[FT_XY, FT_Z] = calcFilterParms(...
    dataProperties.WVL,dataProperties.NA,1,'gauss',...
    dataProperties.sigmaCorrection, ...
    [dataProperties.PIXELSIZE_XY dataProperties.PIXELSIZE_Z]);
FT_XY = FT_XY + fSze_xy/4;
FT_Z = FT_Z + fSze_z/4;

patchXYZ=roundOddOrEven(4*[FT_XY FT_XY FT_Z],'odd','inf');
dataProperties.FILTERPRM = [FT_XY,FT_XY,FT_Z,patchXYZ];

% ===================================================================


spots = spotfind_new(fImg, dataProperties);
spots.start = startTime;
spots.end = clock;

