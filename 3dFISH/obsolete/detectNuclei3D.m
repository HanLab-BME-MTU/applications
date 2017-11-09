function [nuclei] = detectNuclei3D(movieName, Cha, cordVar)

% Treat one nucleus as a single spot and find it
% 2015 Ning

% Take movie data and process to get 3D matrix
Path = strcat('', movieName);
MD = MovieData.load(Path);
xSize = MD.imSize_(1);
ySize = MD.imSize_(2);
nDepth = MD.zSize_;
frameN = 1;

zStack3D = zeros(xSize, ySize, nDepth);
for ii = 1:nDepth
    zStack3D(:,:,ii) = MD.getChannel(Cha).loadImage(frameN,ii);
end
zStack3D=100*zStack3D/max(zStack3D(:));


% Pixel size in um
% XY pixel size 10X calibration error for nikon widefield
dataProperties.PIXELSIZE_XY = MD.pixelSize_/1000*10;
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
gausspar = [fSze_xy/4,fSze_xy/4,fSze_z/4,fSze_xy,fSze_xy,30];
% fImg = filtermovie(zStack3D,gausspar,0);


spotsNum = size(cordVar.sp,2);
nuclei(spotsNum).center = [];
nuclei(spotsNum).cube = [];
for i = 1:spotsNum
    nuclei(i).center = cordVar.sp(i).cord;    
    % Change x y to fit coordinates
    yCtr = round(nuclei(i).center(1));
    xCtr = round(nuclei(i).center(2));
    zCtr = round(nuclei(i).center(3));    
    nuclei(i).cube = stamp3d(zStack3D,gausspar(4:6),[xCtr,yCtr,zCtr],0);
    
%     How to get original nucleus coordinate range??
%     nuclei.nucleus(i).xRange(i).xCord = round((xCtr-fSze_xy/2):(xCtr+fSze_xy/2));
%     nuclei.range(i).yCord = round((yCtr-fSze_xy/2):(yCtr+fSze_xy/2));
%     nuclei.range(i).zCord = round((zCtr-fSze_z/2):(zCtr+fSze_z/2));
end
