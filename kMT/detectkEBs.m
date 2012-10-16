function sisterListEB = detectkEBs(sisterList,sisTrackPairs,tracks,imagesEB,varargin)
%KMTEBSIGNAL measures EB signal at plus-ends of kinetochore microtubules
%
% SYNOPSIS: ebSignalSeries = kmtEBsignal(sisterList,imageEB,varargin)
%
% INPUT sisterList   : Output of groupSisters.
%       sisTrackPairs: Output of groupSisters.
%       tracks       : Output of trackCloseGapsKalmanSparse, in structure
%                      format.
%       imagesEB     : Stack of EB images.
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

%number of frames
imagesEB = double(imagesEB);
[imSize1,imSize2,numFrames] = size(imagesEB);

%% MEASURE EB SIGNAL

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
            
            %output
            sisterListEB(iSis).kEBcoords1 = coords1;
            sisterListEB(iSis).kEBcoords2 = coords2;
            sisterListEB(iSis).kEBamp1 = [signalEB(:,1) NaN(numFrames,2)];
            sisterListEB(iSis).kEBamp2 = [signalEB(:,2) NaN(numFrames,2)];
            
        end
        
    end
    
end

%subtract minimum signal as representative of background
ampMin = vertcat([sisterListEB.kEBamp1; sisterListEB.kEBamp2]);
ampMin = min(ampMin(:))-10;
for iSis = 1 : numSisters
    sisterListEB(iSis).kEBamp1 = sisterListEB(iSis).kEBamp1 - ampMin;
    sisterListEB(iSis).kEBamp2 = sisterListEB(iSis).kEBamp2 - ampMin;
end

%% OUTPUT

