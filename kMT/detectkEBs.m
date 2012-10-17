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
%                           -'radiusEB' [Default: 7].
%
% OUTPUT sisterListEB: Same as input sisterList, but with additional fields
%                      with EB signal information.
%
% created by: Khuloud Jaqaman, October 2012

%% TEST INPUT & READ PARAMETERS

% Input check
assert(isvector(sisterList) && isstruct(sisterList(1)));
assert(ismatrix(sisTrackPairs));
assert(isvector(tracks) && isstruct(tracks(1)));
ip = inputParser;
ip.addParamValue('radiusEB',7,@isscalar);
ip.parse(varargin{:});

% Get parameters
radiusEB = ip.Results.radiusEB;

%number of sisters
numSisters = length(sisterList);
sisTrackPairs = sisTrackPairs(:,1:2);

%remove tracks in sister pairs from list of tracks
numTracks = length(tracks);
indxKeep = setdiff(1:numTracks,sisTrackPairs(:));
tracks = tracks(indxKeep);
numTracks = length(tracks);

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
for iSis = 1 : numSisters
    
    %get coordinates
    coords1 = sisterList(iSis).coords1(:,1:3);
    coords2 = sisterList(iSis).coords2(:,1:3);
    
    %measure EB signal around coordinates in each frame
    signalEB = NaN(numFrames,2);
    for iFrame = 1 : numFrames
        
        %if sisters exist in this frame
        if ~isnan(coords1(iFrame,1))
            
            %make mask around sister 1
            maskEB = zeros(imSize1,imSize2);
            maskEB(round(coords1(iFrame,2)),round(coords1(iFrame,1))) = 1;
            maskEB = imdilate(maskEB,SE);
            maskEB(maskEB==0) = NaN;
            
            %read EB signal
            imageTimesMask = imagesEB(:,:,iFrame) .* maskEB;
            signalEB(iFrame,1) = nanmean(imageTimesMask(:));
            
            %make mask around sister 2
            maskEB = zeros(imSize1,imSize2);
            maskEB(round(coords2(iFrame,2)),round(coords2(iFrame,1))) = 1;
            maskEB = imdilate(maskEB,SE);
            maskEB(maskEB==0) = NaN;
            
            %read EB signal
            imageTimesMask = imagesEB(:,:,iFrame) .* maskEB;
            signalEB(iFrame,2) = nanmean(imageTimesMask(:));
            
        end
        
    end
    
    %subtract background from EB signal
    signalEB = signalEB - bkgSignal;
    
    %if the signal after background subtraction is negative, this means
    %there is no comet
    %indicate that with NaN
    indxNoComet = find(signalEB(:,1)<=0);
    signalEB(indxNoComet,1) = NaN;
    coords1(indxNoComet,:) = NaN;
    indxNoComet = find(signalEB(:,2)<=0);
    signalEB(indxNoComet,2) = NaN;
    coords2(indxNoComet,:) = NaN;
    
    %output
    sisterListEB(iSis).kEBcoords1 = coords1;
    sisterListEB(iSis).kEBcoords2 = coords2;
    sisterListEB(iSis).kEBamp1 = [signalEB(:,1) NaN(numFrames,2)];
    sisterListEB(iSis).kEBamp2 = [signalEB(:,2) NaN(numFrames,2)];
    
end

%% OUTPUT

