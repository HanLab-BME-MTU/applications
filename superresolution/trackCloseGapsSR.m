function [tracksFinal,errFlag] = trackCloseGapsSR(movieInfo,parameterSet,...
    probDim,saveResults,verbose)
%TRACKCLOSEGAPSSR links features between frames and closes gaps for super-resolution applications
%
%SYNOPSIS [tracksFinal,errFlag] = trackCloseGapsSR(movieInfo,costMatrices,...
%    gapCloseParam,probDim,saveResults,verbose)
%
%INPUT  movieInfo    : Array of size equal to the number of frames in a
%                      movie, containing at least the fields:
%             .xCoord      : x-coordinates of detected features. 
%                            1st column: value, 2nd column: standard
%                            deviation (zeros if not available).
%             .yCoord      : y-coordinates of detected features.
%                            1st column: value, 2nd column: standard
%                            deviation (zeros if not available).
%             .zCoord      : z-coordinates of detected features.
%                            1st column: value, 2nd column: standard
%                            deviation (zeros if not available).
%                            Optional. Skipped if problem is 2D. Default: zeros.
%             .amp         : "Intensities" of detected features.
%                            1st column: values (ones if not available),
%                            2nd column: standard deviation (zeros if not
%                            available).
%       parameterSet : Structure with fields:
%             .serachRadius    : Search radius for linking and gap closing.
%             .timeWindow      : Time window for gap closing.
%             .gapPenalty      : Gap penaly when gap closing.
%       probDim      : Problem dimensionality. 2 (for 2D) or 3 (for 3D).
%                      Optional. If not input, dimensionality will be
%                      derived from movieInfo.
%       saveResults  : 0 if no saving is requested.
%                      If saving is requested, structure with fields:
%           .dir          : Directory where results should be saved.
%                           Optional. Default: current directory.
%           .filename     : Name of file where results should be saved.
%                      Or []. Default: trackedFeatures in directory
%                      where run is initiated.
%                      Whole structure optional.
%       verbose      : 1 to show calculation progress, 0 otherwise.
%                      Optional. Default: 1.
%
%       All optional variables can be entered as [] to use default values.
%
%OUTPUT tracksFinal   : Structure array where each element corresponds to a 
%                       compound track. Each element contains the following 
%                       fields:
%           .tracksFeatIndxCG: Connectivity matrix of features between
%                              frames, after gap closing. Number of rows
%                              = number of track segments in compound
%                              track. Number of columns = number of frames
%                              the compound track spans. Zeros indicate
%                              frames where track segments do not exist
%                              (either because those frames are before the
%                              segment starts or after it ends, or because
%                              of losing parts of a segment.
%           .tracksCoordAmpCG: The positions and amplitudes of the tracked
%                              features, after gap closing. Number of rows
%                              = number of track segments in compound
%                              track. Number of columns = 8 * number of
%                              frames the compound track spans. Each row
%                              consists of
%                              [x1 y1 z1 a1 dx1 dy1 dz1 da1 x2 y2 z2 a2 dx2 dy2 dz2 da2 ...]
%                              NaN indicates frames where track segments do
%                              not exist, like the zeros above.
%           .seqOfEvents     : Matrix with number of rows equal to number
%                              of events happening in a compound track and 4
%                              columns:
%                              1st: Frame where event happens;
%                              2nd: 1 = start of track segment, 2 = end of track segment;
%                              3rd: Index of track segment that ends or starts;
%                              4th: NaN = start is a birth and end is a death,
%                                   number = start is due to a split, end
%                                   is due to a merge, number is the index
%                                   of track segment for the merge/split.
%       errFlag       : 0 if function executes normally, 1 otherwise.
%
%Khuloud Jaqaman, August 2011

%% Output

tracksFinal    = [];
errFlag        =  0;

%% Input

%check whether correct number of input arguments was used
if nargin < 2
    disp('--trackCloseGapsSR: Incorrect number of input arguments!');
    errFlag = 1;
    return
end

%get number of frames in movie
numFrames = length(movieInfo);

%check whether z-coordinates were input, making problem potentially 3D
if isfield(movieInfo,'zCoord')
    probDimT = 3;
else
    probDimT = 2;
end

%assign problem dimensionality if not input
if nargin < 3 || isempty(probDim)
    probDim = probDimT;
else
    if probDim == 3 && probDimT == 2
        disp('--trackCloseGapsSR: Inconsistency in input. Problem 3D but no z-coordinates.');
        errFlag = 1;
    end
end

%determine where to save results
if nargin < 4 || isempty(saveResults) %if nothing was input
    saveResDir = pwd;
    saveResFile = 'trackedFeatures';
    saveResults.dir = pwd;
else
    if isstruct(saveResults)
        if ~isfield(saveResults,'dir') || isempty(saveResults.dir)
            saveResDir = pwd;
        else
            saveResDir = saveResults.dir;
        end
        if ~isfield(saveResults,'filename') || isempty(saveResults.filename)
            saveResFile = 'trackedFeatures';
        else
            saveResFile = saveResults.filename;
        end
    else
        saveResults = 0;
    end
end

%check whether verbose
if nargin < 5 || isempty(verbose)
    verbose = 1;
end

%exit if there are problems with input
if errFlag
    disp('--trackCloseGapsSR: Please fix input parameters.');
    return
end

%% preamble

%get tracking parameters
searchRadius = parameterSet.searchRadius;
timeWindow = parameterSet.timeWindow;
gapPenalty = parameterSet.gapPenalty;

%make sure that timeWindow is not equal to 0
%set to 1 in this case, in order to not have any gap closing
if timeWindow == 0
    timeWindow = 1;
end

%get number of features in each frame
if ~isfield(movieInfo,'num')
    for iFrame = 1 : numFrames
        movieInfo(iFrame).num = size(movieInfo(iFrame).xCoord,1);
    end
end

%collect coordinates and their std in one matrix in each frame
if ~isfield(movieInfo,'allCoord')
    switch probDim
        case 2
            for iFrame = 1 : numFrames
                movieInfo(iFrame).allCoord = [movieInfo(iFrame).xCoord ...
                    movieInfo(iFrame).yCoord];
            end
        case 3
            for iFrame = 1 : numFrames
                movieInfo(iFrame).allCoord = [movieInfo(iFrame).xCoord ...
                    movieInfo(iFrame).yCoord movieInfo(iFrame).zCoord];
            end
    end
end

%remove empty frames in the beginning and the end and keep the information
%for later
emptyStart = 0;
numFeatures = vertcat(movieInfo.num);
emptyFrames = find(numFeatures == 0);
if ~isempty(emptyFrames)
    findEmpty = emptyFrames(1) == 1;
else%                           
    findEmpty = 0;
end
while findEmpty
    emptyStart = emptyStart + 1;
    numFeatures = numFeatures(2:end);
    emptyFrames = find(numFeatures == 0);
    if ~isempty(emptyFrames)
        findEmpty = emptyFrames(1) == 1;
    else
        findEmpty = 0;
    end
end
emptyEnd = 0;
numFeatures = vertcat(movieInfo.num);
emptyFrames = find(numFeatures == 0);
if ~isempty(emptyFrames)
    findEmpty = emptyFrames(end) == length(numFeatures);
else
    findEmpty = 0;
end
while findEmpty
    emptyEnd = emptyEnd + 1;
    numFeatures = numFeatures(1:end-1);
    emptyFrames = find(numFeatures == 0);
    if ~isempty(emptyFrames)
        findEmpty = emptyFrames(end) == length(numFeatures);
    else
        findEmpty = 0;
    end
end
movieInfo = movieInfo(emptyStart+1:numFrames-emptyEnd);
numFramesEff = length(movieInfo);

if numFramesEff == 0
    disp('Empty movie. Nothing to track.');
    return
end

%% Link between frames

%get initial track segments by linking features between consecutive frames
if verbose
    disp('Linking features ...');
end
[tracksFeatIndxLink,tracksCoordAmpLink,errFlag] = linkFeaturesSR(...
    movieInfo,searchRadius,probDim,verbose);

%get track segment start and end times
trackSEL = getTrackSEL(tracksCoordAmpLink);
trackStartTime = trackSEL(:,1);
trackEndTime   = trackSEL(:,2);
clear trackSEL

%get number of tracks
numTracksLink = length(trackStartTime);

%% Close gaps

%if there are gaps to close (i.e. if there are tracks that start after the
%first frame and tracks that end before the last frame) ...
if any(trackStartTime > 1) && any(trackEndTime < numFramesEff)

    if verbose
        fprintf('Closing gaps (%d starts and %d ends) ...\n',...
            length(find(trackStartTime>1)),length(find(trackEndTime<numFramesEff)));
    end

    %initialize progress display
    if verbose
        progressText(0,'Gap closing');
    end

    %calculate the cost matrix, which already includes the
    %costs of birth and death
    [costMat,nonlinkMarker,errFlag] = costMatCloseGapsSR(tracksCoordAmpLink,...
    trackStartTime,trackEndTime,searchRadius,timeWindow,gapPenalty,probDim);

    %if there are possible links ...
    if any(isfinite(nonzeros(costMat)))

        %link tracks based on this cost matrix, allowing for birth and death
        [link12,link21] = lap(costMat,nonlinkMarker);
        link12 = double(link12);
        link21 = double(link21);

        %put the indices of all tracks from linking in one vector
        tracks2Link = (1:numTracksLink)';
        tracksRemaining = tracks2Link;

        %reserve memory space for matrix showing track connectivity
        compoundTrack = zeros(numTracksLink,600);

        %initialize compTrackIndx
        compTrackIndx = 0;

        while ~isempty(tracksRemaining)

            %update compound track index by 1
            compTrackIndx = compTrackIndx + 1;

            %take first track as a seed to build a compound track with
            %closed gaps and merges/splits
            trackSeed = tracksRemaining(1);
            seedLength = 1;
            seedLengthOld = 0; %dummy just to get into the while loop

            %while current seed contains more tracks than previous seed, i.e.
            %whie new track segments are still being added to the compound
            %track
            while seedLength > seedLengthOld

                %store current seed for later comparison
                seedLengthOld = seedLength;

                %find tracks connected to ends of seed tracks
                tmpTracks = link12(trackSeed);
                trackLink2End = tmpTracks(tmpTracks <= numTracksLink); %starts linked to ends
                trackMerge = [];

                %find tracks connected to starts of seed tracks
                tmpTracks = link21(trackSeed);
                trackLink2Start = tmpTracks(tmpTracks <= numTracksLink); %ends linked to starts
                trackSplit = [];

                %put all tracks together as the new seed
                trackSeed = [trackSeed; trackLink2End; trackLink2Start; ...
                    trackMerge; trackSplit];

                %remove repetitions and arrange tracks in ascending order
                trackSeed = unique(trackSeed);

                %get number of tracks in new seed
                seedLength = length(trackSeed);

            end %(while length(trackSeed) > length(trackSeedOld))

            %expand trackSeed to reserve memory for connetivity information
            trackSeedConnect = [trackSeed zeros(seedLength,2)];

            %store the tracks that the ends of the seed tracks are linked to,
            %and indicate whether it's an end-to-start link (+ve) or a merge (-ve)
            tmpTracks = link12(trackSeed);
            tmpTracks(tmpTracks > numTracksLink) = NaN;
            trackSeedConnect(:,2) = tmpTracks;

            %store the tracks that the starts of the seed tracks are linked to,
            %and indicate whether it's a start-to-end link (+ve) or a split (-ve)
            tmpTracks = link21(trackSeed);
            tmpTracks(tmpTracks > numTracksLink) = NaN;
            trackSeedConnect(:,3) = tmpTracks;

            %store tracks making up this compound track and their connectivity
            compoundTrack(compTrackIndx,1:3*seedLength) = reshape(...
                trackSeedConnect,3*seedLength,1)';

            %in the list of all tracks, indicate that these tracks have
            %been taken care of by placing NaN instead of their number
            tracks2Link(trackSeed) = NaN;

            %retain only tracks that have not been linked to anything yet
            tracksRemaining = tracks2Link(~isnan(tracks2Link));

        end %(while ~isempty(tracksRemaining))

        %remove empty rows
        maxValue = max(compoundTrack,[],2);
        compoundTrack = compoundTrack(maxValue > 0,:);

        %determine number of tracks after gap closing
        numTracksCG = size(compoundTrack,1);

        %reserve memory for structure storing tracks after gap closing
        tracksFinal = repmat(struct('tracksFeatIndxCG',[],...
            'tracksCoordAmpCG',[],'seqOfEvents',[]),numTracksCG,1);

        %go over all compound tracks
        for iTrack = 1 : numTracksCG

            %get indices of tracks from linking making up current compound track
            %determine their number and connectivity
            trackSeedConnect = compoundTrack(iTrack,:)';
            trackSeedConnect = trackSeedConnect(trackSeedConnect ~= 0);
            seedLength = length(trackSeedConnect)/3; %number of segments making current track
            trackSeedConnect = reshape(trackSeedConnect,seedLength,3);

            %get their start times
            segmentStartTime = trackStartTime(trackSeedConnect(:,1));

            %arrange segments in ascending order of their start times
            [segmentStartTime,indxOrder] = sort(segmentStartTime);
            trackSeedConnect = trackSeedConnect(indxOrder,:);

            %get the segments' end times
            segmentEndTime = trackEndTime(trackSeedConnect(:,1));

            %calculate the segments' positions in the matrix of coordinates and
            %amplitudes
            segmentStartTime8 = 8 * (segmentStartTime - 1) + 1;
            segmentEndTime8   = 8 * segmentEndTime;

            %instead of having the connectivity in terms of the original track
            %indices, have it in terms of the indices of this subset of tracks
            %(which are arranged in ascending order of their start times)
            for iSeed = 1 : seedLength
                value = trackSeedConnect(iSeed,2);
                if value > 0
                    trackSeedConnect(iSeed,2) = find(trackSeedConnect(:,1) == ...
                        value);
                elseif value < 0
                    trackSeedConnect(iSeed,2) = -find(trackSeedConnect(:,1) == ...
                        -value);
                end
                value = trackSeedConnect(iSeed,3);
                if value > 0
                    trackSeedConnect(iSeed,3) = find(trackSeedConnect(:,1) == ...
                        value);
                elseif value < 0
                    trackSeedConnect(iSeed,3) = -find(trackSeedConnect(:,1) == ...
                        -value);
                end
            end

            %get track information from the matrices storing linking information
            tracksFeatIndxCG = tracksFeatIndxLink(trackSeedConnect(:,1),:);
            tracksCoordAmpCG = tracksCoordAmpLink(trackSeedConnect(:,1),:);
            
            %convert zeros to NaNs where approriate for the case of sparse
            %matrices
            if issparse(tracksCoordAmpCG)
                
                %convert sparse to full
                tracksCoordAmpCG = full(tracksCoordAmpCG);
                
                %go over all the rows in this compound track
                for iRow = 1 : size(tracksCoordAmpCG,1)
                    
                    %find all the zero entries
                    colZero = find(tracksCoordAmpCG(iRow,:)==0);
                    colZero = colZero(:)';
                    
                    %find the columns of the x-coordinates corresponding to
                    %the zero columns
                    xCoordCol = colZero - mod(colZero-1,8*ones(size(colZero)));
                    
                    %keep only the columns whose x-coordinate is zero as
                    %well
                    colZero = colZero(tracksCoordAmpCG(iRow,xCoordCol)==0);
                    
                    %replace zero with NaN in the surviving columns
                    tracksCoordAmpCG(iRow,colZero) = NaN;
                    
                end
                
            end

            %perform all gap closing links and modify connectivity accordingly
            %go over all starts in reverse order
            for iSeed = seedLength : -1 : 2

                %find the track this track might be connected to
                track2Append = trackSeedConnect(iSeed,3);

                %if there is a track (which is not a split)
                if track2Append > 0

                    %put track information in the relevant row
                    tracksFeatIndxCG(track2Append,segmentStartTime(iSeed):...
                        segmentEndTime(iSeed)) = tracksFeatIndxCG(iSeed,...
                        segmentStartTime(iSeed):segmentEndTime(iSeed));
                    tracksFeatIndxCG(iSeed,:) = 0;
                    tracksCoordAmpCG(track2Append,segmentStartTime8(iSeed):...
                        segmentEndTime8(iSeed)) = tracksCoordAmpCG(iSeed,...
                        segmentStartTime8(iSeed):segmentEndTime8(iSeed));
                    tracksCoordAmpCG(iSeed,:) = NaN;

                    %update segment information
                    segmentEndTime(track2Append) = segmentEndTime(iSeed);
                    segmentEndTime8(track2Append) = segmentEndTime8(iSeed);
                    segmentEndTime(iSeed) = NaN;
                    segmentEndTime8(iSeed) = NaN;
                    segmentStartTime(iSeed) = NaN;
                    segmentStartTime8(iSeed) = NaN;

                    %update connectivity
                    trackSeedConnect(track2Append,2) = trackSeedConnect(iSeed,2);
                    trackSeedConnect(trackSeedConnect(:,2) == iSeed,2) = track2Append;
                    trackSeedConnect(trackSeedConnect(:,3) == iSeed,3) = track2Append;
                    trackSeedConnect(trackSeedConnect(:,2) == -iSeed,2) = -track2Append;
                    trackSeedConnect(trackSeedConnect(:,3) == -iSeed,3) = -track2Append;

                end %(if track2Append > 0)

            end %(for iSeed = seedLength : -1 : 2)

            %find rows that are not empty
            maxValue = max(tracksFeatIndxCG,[],2);
            rowsNotEmpty = find(maxValue > 0);

            %remove empty rows
            tracksFeatIndxCG = tracksFeatIndxCG(rowsNotEmpty,:);
            tracksCoordAmpCG = tracksCoordAmpCG(rowsNotEmpty,:);
            segmentEndTime   = segmentEndTime(rowsNotEmpty);
            segmentStartTime = segmentStartTime(rowsNotEmpty);
            trackSeedConnect = trackSeedConnect(rowsNotEmpty,:);

            %update connectivity accordingly
            %by now, only merges and splits are left - thus no need for minus
            %sign to distinguish them from closed gaps
            for iSeed = 1 : length(rowsNotEmpty)
                trackSeedConnect(trackSeedConnect(:,2) == -rowsNotEmpty(...
                    iSeed),2) = iSeed;
                trackSeedConnect(trackSeedConnect(:,3) == -rowsNotEmpty(...
                    iSeed),3) = iSeed;
            end

            %determine new "seedLength"
            seedLength = length(rowsNotEmpty);

            %store the sequence of events of this track
            seqOfEvents = [segmentStartTime ones(seedLength,1) ...
                (1:seedLength)' trackSeedConnect(:,3); ...
                segmentEndTime 2*ones(seedLength,1) ...
                (1:seedLength)' trackSeedConnect(:,2)];

            %sort sequence of events in ascending order of time
            [tmp,indxOrder] = sort(seqOfEvents(:,1));
            seqOfEvents = seqOfEvents(indxOrder,:);

            %add 1 to the times of merges
            indx = find(~isnan(seqOfEvents(:,4)) & seqOfEvents(:,2) == 2);
            seqOfEvents(indx,1) = seqOfEvents(indx,1) + 1;

            %find the frame where the compound track starts and the frames
            %where it ends
            frameStart = seqOfEvents(1,1);
            frameEnd   = seqOfEvents(end,1);

            %store final tracks, removing frames before anything happens and
            %after everything happens
            tracksFinal(iTrack).tracksFeatIndxCG = tracksFeatIndxCG(:,...
                frameStart:frameEnd);
            tracksFinal(iTrack).tracksCoordAmpCG = tracksCoordAmpCG(:,...
                8*(frameStart-1)+1:8*frameEnd);
            tracksFinal(iTrack).seqOfEvents = seqOfEvents;

        end %(for iTrack = 1 : numTracksCG)
        
    else %if there are no possible links

        if verbose
            disp('No gaps to close!');
        end

        %convert matrix of tracks into structure
        tracksFinal = convertMat2Struct(tracksCoordAmpLink,tracksFeatIndxLink);

    end %(if any(~isfinite(nonzeros(costMat))))

    %display elapsed time
    if verbose
        progressText(1,'Gap closing');
    end

else %if there are no gaps to close

    if verbose
        disp('No gaps to close!');
    end

    %convert matrix of tracks into structure
    tracksFinal = convertMat2Struct(tracksCoordAmpLink,tracksFeatIndxLink);

end %(if any(trackStartTime > 1) && any(trackEndTime < numFramesEff)

%shift time if any of the initial frames are empty
for iTrack = 1 : length(tracksFinal)
    tracksFinal(iTrack).seqOfEvents(:,1) = tracksFinal(iTrack).seqOfEvents(:,1) + emptyStart;
end

%% Save results

if isstruct(saveResults)
    save([saveResDir filesep saveResFile],'costMatrices','gapCloseParam',...
        'kalmanFunctions','tracksFinal','kalmanInfoLink');
end


%% Subfunction1

function tracksFinal = convertMat2Struct(tracksCoordAmpLink,tracksFeatIndxLink)

%get number of tracks
numTracks = size(tracksCoordAmpLink,1);

%reserve memory for structure storing tracks
tracksFinal = repmat(struct('tracksFeatIndxCG',[],...
    'tracksCoordAmpCG',[],'seqOfEvents',[]),numTracks,1);

%get the start and end time of tracks
trackSEL = getTrackSEL(tracksCoordAmpLink);

%go over all tracks and store information
for iTrack = 1 : numTracks
   
    %track start time and end time
    startTime = trackSEL(iTrack,1);
    endTime = trackSEL(iTrack,2);
    
    %feature indices
    tracksFinal(iTrack).tracksFeatIndxCG = tracksFeatIndxLink(iTrack,startTime:endTime);
    
    %feature coordinates and amplitudes
    tracksFinal(iTrack).tracksCoordAmpCG = full(tracksCoordAmpLink(iTrack,...
        (startTime-1)*8+1:endTime*8));
    
    %sequence of events
    tracksFinal(iTrack).seqOfEvents = [startTime 1 1 NaN; endTime 2 1 NaN];
    
end


%% %%%%% ~~ the end ~~ %%%%%

