function tracksInWindow = assignTracks2Windows(tracksFinal,winPositions,...
    winFrames,assignSegments)
%ASSIGNTRACKS2WINDOWS groups tracks into spatial and temporal windows derived from the cell edge
%
%SYNOPSIS tracksInWindow = assignTracks2Windows(tracksFinal,winPositions,...
%    winFrames,assignSegments)
%
%INPUT  tracksFinal    : The tracks, either in structure format (e.g.
%                        output of trackCloseGapsKalman) or in matrix
%                        format (e.g. output of trackWithGapClosing).
%       winPositions   : The window edges for all time points, as output by
%                        Hunter's old windowing function.
%       winFrames      : The frames at which there are windows.
%       assignSegments : Relevant only for tracks in structure format.
%                        1 to assign track segments, 0 to assign compound
%                        tracks as is.
%                        Optional. Default: 0.
%
%OUTPUT tracksInWindow : Cell array of dimensions (number of bands) x
%                        (number of windows) x (number of window frames-1)
%                        storing the track indices that fall in each window
%                        in each frame.
%
%REMARKS This code is designed for experiments where the particle
%        trajectories are sampled much more frequently than the cell edge.
%        It assumes that particle lifetimes are much shorter than the time
%        between cell edge frames.
%
%        For a different scenario where particle lifetimes are longer than
%        the time between cell edge frames, the tracks should not be
%        grouped like this. Rather, each track should get divided into
%        several segments corresponding to the times between cell edge
%        frames and each track segment should be analyzed separately.
%        Something like that.
%
%Khuloud Jaqaman, May 2010

%% Input

if nargin < 3
    disp('--assignTracks2Windows: Incorrect number of input arguments!');
    return
end

if nargin < 4 || isempty(assignSegments)
    assignSegments = 0;
end
if assignSegments == 1 && ~isstruct(tracksFinal)
    assignSegments = 0;
end

%get number of tracks
numTracksCompound = size(tracksFinal,1);

%get number of windows along the edge and perpendicular to it, and number
%of frames that have windows
[numWinPerp,numWinPara,numWinFrames] = size(winPositions);

%% Pre-processing

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

%% Track assignment into windows

%initialize cell array storing the grouping of tracks into windows for
%each frame range
tracksInWindow = cell(numWinPerp,numWinPara,numWinFrames-1);

%go over all window frames
for iWinFrame = 1 : numWinFrames - 1
    
    %get current frame number and next frame number
    minFrame = winFrames(iWinFrame);
    maxFrame = winFrames(iWinFrame + 1);
    
    %find tracks whose "average" time is in this frame range
    indxFrameRange = find(trackTimeMean>=minFrame & trackTimeMean<maxFrame);
    
    %get the mean positions of these tracks
    xCoordMeanFR = xCoordMean(indxFrameRange);
    yCoordMeanFR = yCoordMean(indxFrameRange);
    
    %go over the windows in this frame
    for iPara = 1 : numWinPara
        for iPerp = 1 : numWinPerp
            
            %if this window has a finite size
            if ~isempty(winPositions(iPerp,iPara,iWinFrame).outerBorder) ...
                    && ~isempty(winPositions(iPerp,iPara,iWinFrame).innerBorder)
                
                %get the window boundaries
                winX = [winPositions(iPerp,iPara,iWinFrame).outerBorder(1,:) ...
                    winPositions(iPerp,iPara,iWinFrame).innerBorder(1,end:-1:1)]';
                winY = [winPositions(iPerp,iPara,iWinFrame).outerBorder(2,:) ...
                    winPositions(iPerp,iPara,iWinFrame).innerBorder(2,end:-1:1)]';
                
                %find the tracks whose "average" position lies in this window
                indxWin = inpolygon(xCoordMeanFR,yCoordMeanFR,winX,winY);
                
                %map back to original track indices
                indxWin = indxFrameRange(indxWin);
                
                %store track indices in cell array
                tracksInWindow{iPerp,iPara,iWinFrame} = indxWin;
                
            else %if this window is collapsed, then there are no tracks in it
                
                tracksInWindow{iPerp,iPara,iWinFrame} = [];
                
            end
            
        end
    end
    
end
