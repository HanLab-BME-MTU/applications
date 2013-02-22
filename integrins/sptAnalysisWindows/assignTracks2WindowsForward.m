function [windowTrackAssign,trackWindowAssign,trackWindowAssignComp,...
    windowTrackAssignExt] = assignTracks2Windows(tracksFinal,winPositions,...
    winFrames,assignSegments)
%ASSIGNTRACKS2WINDOWS groups tracks into spatial and temporal windows derived from the cell edge
%
%SYNOPSIS [windowTrackAssign,trackWindowAssign,trackWindowAssignComp,...
%    windowTrackAssignExt] = assignTracks2Windows(tracksFinal,winPositions,...
%    winFrames,assignSegments)
%
%INPUT  tracksFinal    : The tracks, either in structure format (e.g.
%                        output of trackCloseGapsKalman) or in matrix
%                        format (e.g. output of trackWithGapClosing).
%       winPositions   : A 2D array of the window edges. 
%                        Number of rows = number of window frames. 
%                        Number of columns = number of windows parallel to
%                        Each entry is the output of Hunter's new
%                        windowing software.
%                        Basically, to make this variable, one puts
%                        together the windows of each frame coming out of
%                        the windowing software.
%       winFrames      : The frames at which there are windows.
%       assignSegments : Relevant only for tracks in structure format.
%                        1 to assign track segments, 0 to assign compound
%                        tracks as is.
%                        Optional. Default: 0.
%
%OUTPUT windowTrackAssign: 3D Cell array of dimensions (number of bands) x
%                          (number of slices) x (number of window frames-1)
%                          storing for each window in each frame the track
%                          indices that fall in it.
%       trackWindowAssign: (Number of tracks/track segments) x 3 array
%                          storing for each track/track segment the window
%                          it falls in, indicated by band number, slice
%                          number and frame. trackWindowAssign and
%                          windowTrackAssign are essentially the inverse of
%                          each other.
%       trackWindowAssignComp: Structure array with field "assignment"
%                          storing for each compound track its window
%                          assignment, equivalent to trackWindowAssign. A
%                          "vertcat" of the field "assignment" gives an
%                          array = trackWindowAssign.
%       windowTrackAssignExt: 4D Cell array of dimensions (number of bands) x
%                          (number of slices) x (number of window frames-1)
%                          x (number of window frames-1) storing for each
%                          window in each frame the track indices that fall
%                          in it, not only in its proper frame range but
%                          throughout all frame ranges, as indicated by
%                          the 4th dimension.
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

%generate winFrameMin and winFramesMax for parfor loop
winFramesMin = winFrames(1:end-1);
winFramesMax = winFrames(2:end);

%% Pre-processing

%get number of tracks
numTracksCompound = size(tracksFinal,1);

%get number of frames that have windows and number of windows parallel to
%the edge
[numWinFrames,numWinPara] = size(winPositions);

%find number of windows perpendicular to the edge
nBands = cellfun(@(x)(numel(x)),winPositions);
numWinPerp = max(nBands(:));

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

%initialize cell arrays storing the grouping of tracks into windows for
%each window frame range
windowTrackAssign = cell(numWinPerp,numWinPara,numWinFrames-1);
windowTrackAssignExt = cell(numWinPerp,numWinPara,numWinFrames-1,numWinFrames-1);

%also initialize numTracks x 3 array storing for each track its window
%assignment
trackWindowAssign = NaN(numTracks,3);

%group tracks based on what window frame range they fall into
trackGroup = repmat(struct('indx',[]),numWinFrames-1,1);
parfor iWinFrame = 1 : numWinFrames - 1
   
    %get current spt frame number and next spt frame number
    minFrame = winFramesMin(iWinFrame);
    maxFrame = winFramesMax(iWinFrame);
        
    %find tracks whose "average" time is between minFrame and maxFrame
    indxFrameRange = find(trackTimeMean>=minFrame & trackTimeMean<maxFrame);

    %store this information for later use
    trackGroup(iWinFrame).indx = indxFrameRange;
    
end

%assign tracks to windows
%this is the extended assignement over all frames ranges
parfor iWinFrameExt = 1 : numWinFrames - 1 %this loops over window frames to fetch tracks
    
    %find tracks in this time window
    indxFrameRange = trackGroup(iWinFrameExt).indx;
    
    %get the mean positions of these tracks
    xCoordMeanFR = xCoordMean(indxFrameRange); %#ok<PFBNS>
    yCoordMeanFR = yCoordMean(indxFrameRange); %#ok<PFBNS>
    
    assignMatTmp = cell(numWinPerp,numWinPara,numWinFrames-1);
    for iWinFrame = 1 : numWinFrames - 1 %this loops over window frames to fetch windows
        
        %go over windows in this frame
        for iPara = 1 : numWinPara %slices
            for iPerp = 1 : nBands(iWinFrame,iPara) %#ok<PFBNS> %bands
                
                %if the window in this frame has a finite size
                if ~isempty(winPositions{iWinFrame,iPara}{iPerp}) %#ok<PFBNS>
                    
                    %get the window boundaries
                    windowsPoly = [winPositions{iWinFrame,iPara}{iPerp}{:}];
                    winX = windowsPoly(1,:);
                    winY = windowsPoly(2,:);
                    
                    %find the tracks whose "average" position lies in this
                    %window
                    indxWin = inpolygon(xCoordMeanFR,yCoordMeanFR,winX,winY);
                    
                    %map back to original track indices
                    indxWin = indxFrameRange(indxWin);
                    
                else %if window is collapsed, then there are no tracks in it
                    indxWin = [];
                end
                
                %store extended assignment in temporary matrix
                assignMatTmp{iPerp,iPara,iWinFrame} = indxWin;   
                
            end
        end
        
    end %(for iWinFrame = 1 : numWinFrames - 1)
    
    %store extended assignment in output matrix
    windowTrackAssignExt(:,:,:,iWinFrameExt) = assignMatTmp;
    
end %(for iWinFrameExt = 1 : numWinFrames - 1)

%assign tracks to windows in their proper frame range
for iWinFrame = 1 : numWinFrames - 1
    
    %store track indices in window cell array
    windowTrackAssign(:,:,iWinFrame) = windowTrackAssignExt(:,:,iWinFrame,iWinFrame);
    
    %store window indices for each track
    for iPara = 1 : numWinPara %slices
        for iPerp = 1 : nBands(iWinFrame,iPara) %bands
            indxWin = windowTrackAssign{iPerp,iPara,iWinFrame};
            trackWindowAssign(indxWin,:) = repmat([iPerp iPara iWinFrame],length(indxWin),1);
        end
    end
    
end

%store results in form of compound tracks structure
trackWindowAssignComp = repmat(struct('assignment',[]),numTracksCompound,1);
if assignSegments
    iSeg = 0;
    for iTrack = 1 : numTracksCompound
        numSeg = size(tracksFinal(iTrack).tracksCoordAmpCG,1);
        trackWindowAssignComp(iTrack).assignment = trackWindowAssign(iSeg+1:iSeg+numSeg,:);
        iSeg = iSeg + numSeg;
    end
else
    for iTrack = 1 : numTracks
        trackWindowAssignComp(iTrack).assignment = trackWindowAssign(iTrack,:);
    end
end

