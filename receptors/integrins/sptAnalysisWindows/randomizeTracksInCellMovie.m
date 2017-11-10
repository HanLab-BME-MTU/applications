function tracksFinal = randomizeTracksInCellMovie(tracksFinal,firstMaskFile,...
    maskFrames)
%RANDOMIZETRACKSINCELLMOVIE randomizes the spatial and temporal position of tracks in a movie within the cell
%
%SYNOPSIS tracksFinal = randomizeTracksInCellMovie(tracksIn,firstMaskFile,...
%    maskFrames)
%
%INPUT  tracksFinal    : The original tracks inside the cell mask,
%                        in structure format (i.e. output of
%                        trackCloseGapsKalman).
%       firstMaskFile  : Name of first mask file, including full path.
%       maskFrames     : The frames at which there are masks.
%
%OUTPUT tracksFinal    : The randomized tracks.
%
%REMARKS This code is designed for experiments where the particle
%        trajectories are sampled much more frequently than the cell edge.
%        It assumes that particle lifetimes are much shorter than the time
%        between cell edge frames.
%
%Khuloud Jaqaman, June 2012

%% Input

if nargin < 3
    disp('--randomizeTracksInCellMovie: Incorrect number of input arguments!');
    return
end

%get names of mask files
outFileList = getFileStackNames(firstMaskFile);
numMaskFrames = length(outFileList);

%get image size
mask1 = imread(outFileList{1});
[imSizeY,imSizeX] = size(mask1);
imSizeTot = imSizeX*imSizeY;

%% Randomize time

%get number of tracks
numTracks = size(tracksFinal,1);

%get track start and end times
trackSEL = getTrackSEL(tracksFinal);

%calculate the "average" time at which a track exists
trackTimeMean0 = mean(trackSEL(:,1:2),2);

%assign tracks a random "average" time
rng('shuffle','twister');
trackTimeMean = randi([floor(min(trackTimeMean0)) ceil(max(trackTimeMean0))],numTracks,1) + ...
    (1-mod(trackSEL(:,3),2))/2;
trackTimeShift = trackTimeMean - trackTimeMean0;
for iTrack = 1 : numTracks
    tracksFinal(iTrack).seqOfEvents(:,1) = tracksFinal(iTrack).seqOfEvents(:,1) ...
        + trackTimeShift(iTrack);
end

%% Randomize position

%get the average position of each track
[xCoordMean0,yCoordMean0,xCoordMean,yCoordMean] = deal(NaN(numTracks,1));
for iTrack = 1 : numTracks
    xCoordAll = tracksFinal(iTrack).tracksCoordAmpCG(:,1:8:end);
    yCoordAll = tracksFinal(iTrack).tracksCoordAmpCG(:,2:8:end);
    xCoordMean0(iTrack) = nanmean(xCoordAll(:));
    yCoordMean0(iTrack) = nanmean(yCoordAll(:));
end

%go over all cell masks
for iMaskFrame = 1 : numMaskFrames - 1
    
    %read current mask and next mask
    mask1 = imread(outFileList{iMaskFrame});
    mask2 = imread(outFileList{iMaskFrame+1});
    mask = mask1 | mask2;
    
    %determine how much of the image is covered by the combined mask
    imFrac = length(find(mask)) / imSizeTot;

    %get current frame number and next frame number
    minFrame = maskFrames(iMaskFrame);
    maxFrame = maskFrames(iMaskFrame + 1);
    
    %find tracks whose "average" time is in this frame range
    indxFrameRange = find(trackTimeMean>=minFrame & trackTimeMean<maxFrame);
    numTracksRange = length(indxFrameRange);
    
    %estimate number of random positions to pick
    %double that number in order to avoid looping
    numPos = ceil(2 * numTracksRange / imFrac);
    
    %get random positions
    xMeanTmp = rand(numPos,1)*imSizeX;
    yMeanTmp = rand(numPos,1)*imSizeY;
    
    %keep only those positions inside the cell mask
    linIndx = sub2ind([imSizeY imSizeX],ceil(yMeanTmp),ceil(xMeanTmp));
    indxKeep = find(mask(linIndx));
    xMeanTmp = xMeanTmp(indxKeep);
    yMeanTmp = yMeanTmp(indxKeep);
    
    %loop if necessary until the required number of positions is attained
    while length(xMeanTmp) < numTracksRange
        
        %get additional random positions
        xMeanTmp2 = rand(numPos,1)*imSizeX;
        yMeanTmp2 = rand(numPos,1)*imSizeY;
        
        %keep only those positions inside the cell mask
        linIndx = sub2ind([imSizeY imSizeX],ceil(yMeanTmp2),ceil(xMeanTmp2));
        indxKeep = find(mask(linIndx));
        xMeanTmp = [xMeanTmp; xMeanTmp2(indxKeep)];
        yMeanTmp = [yMeanTmp; yMeanTmp2(indxKeep)];
        
    end
    
    %assign new random average positions to tracks in this frame range
    xCoordMean(indxFrameRange) = xMeanTmp(1:numTracksRange);
    yCoordMean(indxFrameRange) = yMeanTmp(1:numTracksRange);
    
end

%calculate the position shift for each track
xCoordShift = xCoordMean - xCoordMean0;
yCoordShift = yCoordMean - yCoordMean0;

%thus shift the actual coordinates for each track
for iTrack = 1 : numTracks
    tracksFinal(iTrack).tracksCoordAmpCG(:,1:8:end) = ...
        tracksFinal(iTrack).tracksCoordAmpCG(:,1:8:end) + xCoordShift(iTrack);
    tracksFinal(iTrack).tracksCoordAmpCG(:,2:8:end) = ...
        tracksFinal(iTrack).tracksCoordAmpCG(:,2:8:end) + yCoordShift(iTrack);
end

%% ~~~ the end ~~~

