function spots = spotfindNormchannels(movieName)
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

% 07/2014 Ning

% Take movie data and process to get 3D matrix
Path = strcat('/home2/nzhang/', movieName);
MD = MovieData.load(Path);
xSize = MD.imSize_(1);
ySize = MD.imSize_(2);
nDepth = MD.zSize_;
frameN = 1;

% Generate dataProperties.FILTERPRM
dataProperties.NA=1.1000;
dataProperties.WVL = {0.610, 0.509, 0.461};
% Emission wave length of 610nm from mCherry (rfp), 509nm for gfp, 461nm for dapi 
dataProperties.sigmaCorrection=[1 1];
dataProperties.PIXELSIZE_XY = MD.pixelSize_/1000;
dataProperties.PIXELSIZE_Z = MD.pixelSizeZ_/1000;
avgStack3D = zeros(xSize, ySize, nDepth);

for Cha = 1:3    

    [FT_XY, FT_Z] = calcFilterParms(...
        dataProperties.WVL{Cha},dataProperties.NA,1.334,'gauss',...
        dataProperties.sigmaCorrection, ...
        [dataProperties.PIXELSIZE_XY dataProperties.PIXELSIZE_Z]);
    patchXYZ=roundOddOrEven(4*[FT_XY FT_XY FT_Z],'odd','inf');
    dataProperties.FILTERPRM = [FT_XY,FT_XY,FT_Z,patchXYZ];

    zStack3D = zeros(xSize, ySize, nDepth);
    for ii = 1:nDepth
        zStack3D(:,:,ii) = MD.getChannel(Cha).loadImage(frameN,ii);
    end
    zStack3D=100*zStack3D/max(zStack3D(:));
    avgStack3D = (avgStack3D*(Cha-1)+zStack3D)/Cha;
end


% Gaussian filter sigma size
gausspar = ones(1,6)*33;
fImg = filtermovie(avgStack3D,gausspar,0);
spots = spotfind_new(fImg, dataProperties);

% spots = spotfind_new(avgStack3D, dataProperties);
