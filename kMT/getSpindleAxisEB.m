function [spindleAxisVec,poleInfo] = getSpindleAxisEB(imageEB,varargin)
%GETSPINDLEAXIS approximates the spindle axis from EB images in 2D
%
%SYNPOSIS [spindleAxisVec,polePos] = getSpindleAxisEB(imageEB,doPlot)
%
%INPUT imageEB : EB image to be used for spindle axis estimation
%      varargin: Input parameters in param-value pairs:
%           'doPlot' - 0/1 (default: 0).
%
%OUTPUT spindleAxisVec: Array with normalized spindle axis vector per
%                       frame.
%       poleInfo      : Pole positions and amplitudes in each frame, in the
%                       form of movieInfo.
%
%Khuloud Jaqaman, October 2012

%% Input

%input check
ip = inputParser;
ip.addOptional('doPlot',0,@isscalar);
ip.parse(varargin{:});

%get optional parameters
doPlot=ip.Results.doPlot;

%get number of frames
numFrames = size(imageEB,3);

%make sure image is in double format
imageEB = double(imageEB);

%% Determine spindle axis

%average all images
imageEB = mean(imageEB,3);

%find bright spindle pole
prctVal = 100 * (1 - 100/numel(imageEB));
threshPole = prctile(imageEB(:),prctVal);
maskPole = imageEB >= threshPole;

%get the location of the bright spindle pole as the centroid of its mask
%in case there are more than one connected region in the mask, take the
%largest one
stats = regionprops(maskPole,imageEB,'Area','Centroid','MeanIntensity');
maskArea = vertcat(stats.Area);
maskCentroid = vertcat(stats.Centroid);
maskIntensity = vertcat(stats.MeanIntensity);
indxPole1 = find(maskArea==max(maskArea));
pole1Pos = round(maskCentroid(indxPole1,:)); %image coordinates
pole1Amp = maskIntensity(indxPole1);
maskPole(:) = 0;
maskPole(pole1Pos(2),pole1Pos(1)) = 1;

%calculate the distance of each pixel from the bright spindle pole
distFromPole = bwdist(maskPole);

%determine mask for whole spindle
threshSpindle = prctile(imageEB(distFromPole<50),40);
maskSpindle = imageEB >= threshSpindle;
maskSpindle = imfill(maskSpindle,'holes');

%if there is more than one connected region, keep the largest
stats = regionprops(maskSpindle,'Area','PixelIdxList');
maskArea = vertcat(stats.Area);
if length(maskArea) > 1
    maskSpindle(:) = 0;
    maskSpindle(stats(maskArea==max(maskArea)).PixelIdxList) = 1;
end

%keep only distances inside the spindle mask
distFromPole = distFromPole .* maskSpindle;

%filter the distances image to remove any funky structures
distFromPole = filterGauss2D(distFromPole,5);

%find the pixel with the largest distance, this is the location of the
%other pole
[pole2PosY,pole2PosX] = find(distFromPole==max(distFromPole(:))); %image coordinates
pole2Pos = [pole2PosX pole2PosY];
imageCrop = imageEB(pole2PosY-5:pole2PosY+5,pole2PosX-5:pole2PosX+5);
pole2Amp = mean(imageCrop(:));

%calculate the spindle axis vector
spindleAxisVec = pole2Pos - pole1Pos;
spindleAxisVec = spindleAxisVec / norm(spindleAxisVec);

%estimate the average background intensity and subtract it from pole
%amplitudes
SE = strel('square',5);
mask2 = imdilate(maskSpindle,SE);
mask2 = double(mask2) - double(maskSpindle);
mask2(mask2==0) = NaN;
tmp = mask2 .* imageEB;
bkgSignal = nanmean(tmp(:));

%output
spindleAxisVec = repmat(spindleAxisVec,numFrames,1);
polePosBoth = [pole1Pos; pole2Pos];
poleAmpBoth = [pole1Amp; pole2Amp] - bkgSignal;
poleInfo = struct('xCoord',[polePosBoth(:,1) zeros(2,1)],'yCoord',[polePosBoth(:,2) zeros(2,1)],'amp',[poleAmpBoth zeros(2,1)]);
poleInfo = repmat(poleInfo,numFrames,1);

%% Plot

if doPlot
    
    imageEBNorm = (imageEB - min(imageEB(:))) / (max(imageEB(:)) - min(imageEB(:)));
    maskPerim = maskSpindle - imerode(maskSpindle,strel('square',3));
    figure
    imshow(cat(3,maskPerim,imageEBNorm,imageEBNorm),[]);
    hold on
    plot(pole1Pos(1),pole1Pos(2),'rx','MarkerSize',20)
    plot(pole2Pos(1),pole2Pos(2),'rx','MarkerSize',20)
    
end

%% ~~~ the end ~~~