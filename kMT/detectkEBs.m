function sisterListEB = detectkEBs(sisterList,sisTrackPairs,tracks,imagesEB,varargin)
%KMTEBSIGNAL measures EB signal at plus-ends of kinetochore microtubules
%
% SYNOPSIS: sisterListEB = detectkEBs(sisterList,sisTrackPairs,tracks,imagesEB,varargin)
%
% INPUT sisterList   : Output of groupSisters.
%       sisTrackPairs: Output of groupSisters.
%       tracks       : Output of trackCloseGapsKalmanSparse, in structure
%                      format.
%       imagesEB     : Stack of EB images.
%       varargin     : Optional input variables in the form of variable
%                      name/value pairs. Currently this includes:
%                           -'radiusEB' [Default: 5].
%                           -'lengthAlongMT' [Default: 10].
%
% OUTPUT sisterListEB: Same as input sisterList, but with additional fields
%                      with EB signal information.
%
% REMARKS Code currently assumes that spindle poles have been detected and
% thus the pole-kinetochore vectors are known. Must be changed to handle
% case when this is not true. For example, use either only a round area to
% get EB intensity, or make an ellipitical area using the sister-sister
% vector.
%
% created by: Khuloud Jaqaman, October 2012

%% TEST INPUT & READ PARAMETERS

% Input check
assert(isvector(sisterList) && isstruct(sisterList(1)));
assert(ismatrix(sisTrackPairs));
assert(isvector(tracks) && isstruct(tracks(1)));
ip = inputParser;
ip.addParamValue('radiusEB',5,@isscalar);
ip.addParamValue('lengthAlongMT',10,@isscalar);
ip.parse(varargin{:});

% Get parameters
radiusEB = ip.Results.radiusEB;
lengthAlongMT = ip.Results.lengthAlongMT;

%number of sisters
numSister = length(sisterList);
sisTrackPairs = sisTrackPairs(:,1:2);

%remove tracks in sister pairs from list of tracks
numTrack = length(tracks);
indxKeep = setdiff(1:numTrack,sisTrackPairs(:));
tracks = tracks(indxKeep);
numTrack = length(tracks);

%image size and number of frames
imagesEB = double(imagesEB);
[imSize1,imSize2,numFrames] = size(imagesEB);

%% ESTIMATE EB BACKGROUND SIGNAL - SAME AS IN getSpindleAxisEB

%get average EB image
meanImageEB = mean(imagesEB,3);

%find bright spindle pole
prctVal = 100 * (1 - 100/numel(meanImageEB));
threshPole = prctile(meanImageEB(:),prctVal);
maskPole = meanImageEB >= threshPole;

%get the location of the bright spindle pole as the centroid of its mask
%in case there are more than one connected region in the mask, take the
%largest one
stats = regionprops(maskPole,'Area','Centroid');
maskArea = vertcat(stats.Area);
maskCentroid = vertcat(stats.Centroid);
pole1Pos = round(maskCentroid(maskArea==max(maskArea),:)); %image coordinates
maskPole(:) = 0;
maskPole(pole1Pos(2),pole1Pos(1)) = 1;

%calculate the distance of each pixel from the bright spindle pole
distFromPole = bwdist(maskPole);

%determine mask for whole spindle
threshSpindle = prctile(meanImageEB(distFromPole<50),40);
maskSpindle = meanImageEB >= threshSpindle;
maskSpindle = imfill(maskSpindle,'holes');

%if there is more than one connected region, keep the largest
stats = regionprops(maskSpindle,'Area','PixelIdxList');
maskArea = vertcat(stats.Area);
if length(maskArea) > 1
    maskSpindle(:) = 0;
    maskSpindle(stats(maskArea==max(maskArea)).PixelIdxList) = 1;
end

%estimate the average background intensity
SE = strel('square',5);
mask2 = imdilate(maskSpindle,SE);
mask2 = double(mask2) - double(maskSpindle);
mask2(mask2==0) = NaN;
tmp = mask2 .* meanImageEB;
bkgSignal = nanmean(tmp(:));

%% MEASURE KINETOCHORE EB SIGNAL

%output
sisterListEB = sisterList;

%design structuring element
SE = strel('disk',radiusEB,0);

%go over all sisters and get the associated EB signal
for iSis = 1 : numSister
    
    %get coordinates and vectors from poles
    coordsSis = cat(3,sisterList(iSis).coords1(:,1:3),sisterList(iSis).coords2(:,1:3));
    vecFromPole = cat(3,sisterList(iSis).vecFromPole1(:,1:3),sisterList(iSis).vecFromPole2(:,1:3));
    
    %measure EB signal around coordinates in each frame
    signalEB = NaN(numFrames,2);
    for iFrame = 1 : numFrames
        
        %if sisters exist in this frame
        if ~isnan(coordsSis(iFrame,1,1))
            
            %             imageTmp = zeros(imSize1,imSize2,3);
            %             imageTmp(:,:,1) = imagesEB(:,:,iFrame)/max(max(imagesEB(:,:,iFrame)));
            
            for iKin = 1 : 2
                
                %make line to generate mask around kinetochore
                coordsKin = round(coordsSis(iFrame,:,iKin));
                vecKin = vecFromPole(iFrame,:,iKin);
                vecKin = vecKin / normList(vecKin);
                coordsKinExt = round(coordsKin - vecKin*lengthAlongMT);
                coordsLine(:,1) = round(linspace(coordsKin(1),coordsKinExt(1),round(2*lengthAlongMT)));
                coordsLine(:,2) = round(linspace(coordsKin(2),coordsKinExt(2),round(2*lengthAlongMT)));
                coordsLineMat = sub2ind([imSize1 imSize2],coordsLine(:,2),coordsLine(:,1));
                
                %make mask around kinetochore
                maskEB = zeros(imSize1,imSize2);
                maskEB(coordsLineMat) = 1;
                maskEB = imdilate(maskEB,SE);
                %                 maskCirc = imdilate(maskEB,strel('square',3))-maskEB;
                %                 imageTmp(:,:,iKin+1) = maskCirc;
                
                %read EB signal
                maskEB(maskEB==0) = NaN;
                imageTimesMask = imagesEB(:,:,iFrame) .* maskEB;
                signalEB(iFrame,iKin) = nanmean(imageTimesMask(:));
                
            end
            
            %             imshow(imageTmp,[])
            
        end
        
    end
    
    %subtract background from EB signal
    signalEB = signalEB - bkgSignal;
    
    %if the signal after background subtraction is negative, this means
    %there is no comet
    %indicate that with 0 forthe intensity and NaN for the coordinates
    for iKin = 1 : 2
        indxNoComet = find(signalEB(:,iKin)<=0);
        signalEB(indxNoComet,iKin) = 0;
        coordsSis(indxNoComet,:,iKin) = NaN;
    end
    
    %output
    sisterListEB(iSis).kEBcoords1 = coordsSis(:,:,1);
    sisterListEB(iSis).kEBcoords2 = coordsSis(:,:,2);
    sisterListEB(iSis).kEBamp1 = [signalEB(:,1) NaN(numFrames,2)];
    sisterListEB(iSis).kEBamp2 = [signalEB(:,2) NaN(numFrames,2)];
    
end

%% OUTPUT

