function [indxTracksInCellMask] = findTracksInCellMask(tracksFinal,...
    firstMaskFile,maskFrames,assignSegments)
%FINDTRACKSINCELLMASK returns the indices of tracks within the boundaries of the cell mask
%
%SYNOPSIS [indxTracksInCellMask] = findTracksInCellMask(tracksFinal,...
%    firstMaskFile,maskFrames,assignSegments)
%
%INPUT  tracksFinal    : The tracks, either in structure format (e.g.
%                        output of trackCloseGapsKalman) or in matrix
%                        format (e.g. output of trackWithGapClosing).
%       firstMaskFile  : Name of first mask file, including full path.
%       maskFrames     : The frames at which there are masks.
%       assignSegments : Relevant only for tracks in structure format.
%                        1 to assign track segments, 0 to assign compound
%                        tracks as is.
%                        Optional. Default: 0.
%
%OUTPUT indxTracksInCellMask: vector with indices of tracks in cell masks.
%
%REMARKS This code is designed for experiments where the particle
%        trajectories are sampled much more frequently than the cell edge.
%        It assumes that particle lifetimes are much shorter than the time
%        between cell edge frames.
%
%Khuloud Jaqaman, June 2012

%% Input

if nargin < 3
    disp('--findTracksInCellMask: Incorrect number of input arguments!');
    return
end

if nargin < 4 || isempty(assignSegments)
    assignSegments = 0;
end
if assignSegments == 1 && ~isstruct(tracksFinal)
    assignSegments = 0;
end

%get names of mask files
outFileList = getFileStackNames(firstMaskFile);
numMaskFrames = length(outFileList);

%% Pre-processing

%get number of tracks
numTracksCompound = size(tracksFinal,1);

%get track/track segment start and end times
trackSEL = getTrackSEL(tracksFinal,assignSegments);

%get number of tracks/track segments
numTracks = size(trackSEL,1);

%calculate the "average" time at which a track exists
%this will be used to assign tracks to time windows
trackTimeMean = mean(trackSEL(:,1:2),2);

%get average track positions
if isstruct(tracksFinal) %if compound tracks
    
    %initilize mean position
    [xCoordMean,yCoordMean] = deal(NaN(numTracks,1));
    
    %get mean position based on compound tracks or track segments
    if assignSegments
        iSeg = 0;
        for iTrack = 1 : numTracksCompound
            xCoordAll = tracksFinal(iTrack).tracksCoordAmpCG(:,1:8:end);
            yCoordAll = tracksFinal(iTrack).tracksCoordAmpCG(:,2:8:end);
            numSeg = size(xCoordAll,1);
            xCoordMean(iSeg+1:iSeg+numSeg) = nanmean(xCoordAll,2);
            yCoordMean(iSeg+1:iSeg+numSeg) = nanmean(yCoordAll,2);
            iSeg = iSeg + numSeg;
        end
    else
        for iTrack = 1 : numTracks
            xCoordAll = tracksFinal(iTrack).tracksCoordAmpCG(:,1:8:end);
            yCoordAll = tracksFinal(iTrack).tracksCoordAmpCG(:,2:8:end);
            xCoordMean(iTrack) = nanmean(xCoordAll(:));
            yCoordMean(iTrack) = nanmean(yCoordAll(:));
        end
    end
    
else %if simple tracks
    
    %extract x- and y-coordinates
    xCoordAll = tracksFinal(:,1:8:end);
    yCoordAll = tracksFinal(:,2:8:end);
    
    %calculate the average coordinates
    xCoordMean = nanmean(xCoordAll,2);
    yCoordMean = nanmean(yCoordAll,2);
    
end

%% Track classification inside and outside of cell masks

%go over all mask frames - 1
indxKeep = repmat(struct('values',[]),numMaskFrames,1);
for iMaskFrame = 1 : numMaskFrames - 1
    
    %read current mask and next mask
    mask1 = imread(outFileList{iMaskFrame});
    mask2 = imread(outFileList{iMaskFrame+1});
    mask = mask1 | mask2;
    
    %get current frame number and next frame number
    minFrame = maskFrames(iMaskFrame);
    maxFrame = maskFrames(iMaskFrame + 1);
    
    %find tracks whose "average" time is in this frame range
    indxFrameRange = find(trackTimeMean>=minFrame & trackTimeMean<maxFrame);
    
    %get the mean positions of these tracks
    xCoordMeanFR = xCoordMean(indxFrameRange);
    yCoordMeanFR = yCoordMean(indxFrameRange);
    
    %convert the x and y coordinates into a linear index
    linIndx = sub2ind(size(mask),ceil(yCoordMeanFR),ceil(xCoordMeanFR));
    
    %find tracks that lie within mask
    indxKeep(iMaskFrame).values = indxFrameRange(mask(linIndx));
        
end

indxTracksInCellMask = vertcat(indxKeep.values);

%% ~~~ the end ~~~

