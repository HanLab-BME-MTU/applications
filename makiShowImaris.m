function makiShowImaris(dataStruct,select)
%MAKISHOWIMARIS shows mammalina kinetochore data in Imaris
%
% SYNOPSIS: imarisHandle = makiShowImaris(dataStruct,select)
%
% INPUT dataStruct: (opt) data structure as described in makiMakeDataStruct
%                   if empty, guiLoad
%		select: (opt) vector of analysis results to plot. If not
%               not inputed or empty, everything available will be plotted.
%               1st entry: tracks; 2nd entry: sisters, 3rd entry: fitted 
%               plane. Enter 1 for result to be plotted, 0 for result not 
%               to be plotted.
%
% OUTPUT ---
%
% REMARKS in plotting tracks, merges and splits cannot be plotted
%
% created with MATLAB ver.: 7.4.0.287 (R2007a) on Windows_NT
%
% created by: jdorn, kjaqaman
% DATE: 29-Jun-2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% input
if nargin < 1
    dataStruct = [];
end
if isempty(dataStruct)
    dataStruct = makiLoadDataFile;
end

if nargin < 2 || isempty(select)
    select = [1 1 1];
elseif length(select) < 2
    select = [select 1 1];
elseif length(select) < 3
    select = [select 1];
end
    
% turn off property reader warning
warningState = warning;
warning off IMARISIMREAD:NOPROPERTYREADER

% reduce amount of typing
dataProperties = dataStruct.dataProperties;
pixelSize = [dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_XY,...
    dataProperties.PIXELSIZE_Z];

% check for croppping
if isempty(dataProperties.crop)
    dataProperties.crop = zeros(2,3);
end
crop = dataProperties.crop(:,1:3);
isCrop = any(crop,1);
crop(1,~isCrop) = 1;
crop(2,~isCrop) = dataProperties.movieSize(find(~isCrop)); %#ok<FNDSB>

% start imaris
imarisApplication = imarisStartNew;

% load raw movie into imaris. We could do filtered movie, but
% for this, we would have to load frame by frame and do all the
% image properties stuff
imarisApplication.FileOpen(...
    fullfile(dataStruct.rawMoviePath,dataStruct.rawMovieName),...
    'reader=''DeltaVision''');

% check image properties: image should begin at -0.5 pix.
% It would be nice to be able to set the pixelSize to -0.5.
% Instead, we have to read the mins and calculate an offset
zeroOffsetX = imarisApplication.mDataSet.mExtendMinX + 0.5*pixelSize(1);
zeroOffsetY = imarisApplication.mDataSet.mExtendMinY + 0.5*pixelSize(2);
zeroOffsetZ = imarisApplication.mDataSet.mExtendMinZ + 0.5*pixelSize(3);
zeroOffset = [zeroOffsetX zeroOffsetY zeroOffsetZ];

%% detected spots (initCoord)

%get coordinates
initCoord = dataStruct.initCoord;

% make spots object: X,Y,Z,T,r
nTimepoints = dataProperties.movieSize(end);
nSpots = cat(1,initCoord.nSpots);
spots = zeros(sum(nSpots),5);
spots(:,5) = pixelSize(1)*2; % radius in micron
goodTimes = find(nSpots);
nSpotSum = [0;cumsum(nSpots)];
for t = goodTimes'
    % calculate positions in microns. Subtract one voxel from
    % the coords: Imaris starts counting at 0!
    % Use initCoord in pixels to avoid correction problems
    spots(nSpotSum(t)+1:nSpotSum(t+1),1:4) = ...
        [(initCoord(t).allCoordPix(:,[2,1,3])-1 + ...
        repmat(crop(1,[2,1,3])-1,nSpots(t),1)).*...
        repmat(pixelSize,nSpots(t),1) + ...
        repmat(zeroOffset,nSpots(t),1),...
        (t-1)*ones(nSpots(t),1)];
end

% make top-level surpass scene
imaSurpassScene = imarisApplication.mFactory.CreateDataContainer();

% fill surpass scene with light and frame and volume
imaLight = imarisApplication.mFactory.CreateLightSource();
imaSurpassScene.AddChild(imaLight);
imaFrame = imarisApplication.mFactory.CreateFrame();
imaSurpassScene.AddChild(imaFrame);
imaVolume = imarisApplication.mFactory.CreateVolume();
imaSurpassScene.AddChild(imaVolume);

% add surpass scene and set view
imarisApplication.mSurpassScene = imaSurpassScene;
imarisApplication.mViewer = 'eViewerSurpass';

% create spots object
imaSpots = imarisApplication.mFactory.CreateSpots;

% set coords
imaSpots.Set(single(spots(:,1:3)),single(spots(:,4)),single(spots(:,5)));

%assign name
imaSpots.mName = ['Spots (avg: ' num2str(round(mean(nSpots))) ' / frame)'];

% add to scene
imaSurpassScene.AddChild(imaSpots);
           
%% tracks

if select(1)

    %make spots plotted earlier invisible
    imaSpots.mVisible = 0;

    %get tracks from dataStruct
    tracksFinal = dataStruct.tracks;

    %find total number of tracks
    numTracks = length(tracksFinal);

    %find track start times, end times and lifetimes
    trackSEL = getTrackSEL(tracksFinal);

    %find gaps in tracks
    gapInfo = findTrackGaps(tracksFinal);

    %create data container for tracks longer than 90% of movie
    imaTrackGroup90to100 = imarisApplication.mFactory.CreateDataContainer;

    %create data container for tracks between 70% and 90% of movie
    imaTrackGroup70to90 = imarisApplication.mFactory.CreateDataContainer;
    imaTrackGroup70to90.mName = 'tracks - length 70-90%';

    %create data container for tracks between 50% and 70% of movie
    imaTrackGroup50to70 = imarisApplication.mFactory.CreateDataContainer;
    imaTrackGroup50to70.mName = 'tracks - length 50-70%';

    %create data container for tracks shorter than 50% of movie
    imaTrackGroup0to50 = imarisApplication.mFactory.CreateDataContainer;
    imaTrackGroup0to50.mName = 'tracks - length 0-50%';

    %initialize to zero number of tracks in each category
    numTracks90to100 = 0;
    numTracks70to90 = 0;
    numTracks50to70 = 0;
    numTracks0to50 = 0;

    %plot unpaired tracks
    for iTrack = 1 : numTracks

        %create track object
        imaTracks = imarisApplication.mFactory.CreateTrack;

        %get spots belonging to this track, where index is per
        %frame
        spotsIndx = [ones(1,trackSEL(iTrack,1)-1) ...
            tracksFinal(iTrack).tracksFeatIndxCG ...
            ones(1,nTimepoints-trackSEL(iTrack,2))]';

        %locate gaps in this track
        gapsInTrack = gapInfo(gapInfo(:,1)==iTrack,:);

        %calculate cumulative index of spots in order to get spot
        %data from the variable "spots"
        spotsIndx = spotsIndx + nSpotSum(1:end-1);

        %get spot coordinates (some are wrong and will be corrected
        %in the next couple of steps) and define spot sizes 
        spotsCoord = spots(spotsIndx,1:3);
        spotSize = pixelSize(1)*2*ones(nTimepoints,1);

        %for frames before track starts, assign position as that at
        %the start. Make spot size 0
        spotsCoord(1:trackSEL(iTrack,1)-1,1) = spotsCoord(trackSEL(iTrack,1),1);
        spotsCoord(1:trackSEL(iTrack,1)-1,2) = spotsCoord(trackSEL(iTrack,1),2);
        spotsCoord(1:trackSEL(iTrack,1)-1,3) = spotsCoord(trackSEL(iTrack,1),3);
        spotSize(1:trackSEL(iTrack,1)-1) = 0;

        %for frames after track ends, assign position as that at
        %the end. Make spot size 0
        spotsCoord(trackSEL(iTrack,2)+1:end,1) = spotsCoord(trackSEL(iTrack,2),1);
        spotsCoord(trackSEL(iTrack,2)+1:end,2) = spotsCoord(trackSEL(iTrack,2),2);
        spotsCoord(trackSEL(iTrack,2)+1:end,3) = spotsCoord(trackSEL(iTrack,2),3);
        spotSize(trackSEL(iTrack,2)+1:nTimepoints) = 0;

        %in frames where there is a gap, use coordinate of last
        %frame where object is detected. Make spot size half that of a
        %detected spot
        for iGap = 1 : size(gapsInTrack,1)
            spotsCoord(gapsInTrack(iGap,3):gapsInTrack(iGap,3)+gapsInTrack(iGap,4)-1,1) = spotsCoord(gapsInTrack(iGap,3)-1,1);
            spotsCoord(gapsInTrack(iGap,3):gapsInTrack(iGap,3)+gapsInTrack(iGap,4)-1,2) = spotsCoord(gapsInTrack(iGap,3)-1,2);
            spotsCoord(gapsInTrack(iGap,3):gapsInTrack(iGap,3)+gapsInTrack(iGap,4)-1,3) = spotsCoord(gapsInTrack(iGap,3)-1,3);
            spotSize(gapsInTrack(iGap,3):gapsInTrack(iGap,3)+gapsInTrack(iGap,4)-1) = pixelSize(1);
        end

        %set spot coordinates in imaris object
        imaSpotsTrack = imarisApplication.mFactory.CreateSpots;
        imaSpotsTrack.Set(single(spotsCoord),...
            single(0:nTimepoints-1),single(spotSize));

        %define track spots
        imaTracks.SetSpots(imaSpotsTrack);

        %define track edges
        imaTracks.SetEdges(single([(0:nTimepoints-2)' (1:nTimepoints-1)']));

        %add track to relevant data container
        if trackSEL(iTrack,3) >= 0.9*nTimepoints
            imaTracks.SetColor(single(1),single(0),single(0),single(0));
            imaTrackGroup90to100.AddChild(imaTracks);
            numTracks90to100 = numTracks90to100 + 1;
        elseif trackSEL(iTrack,3) >= 0.7*nTimepoints
            imaTracks.SetColor(single(0),single(1),single(0),single(0));
            imaTrackGroup70to90.AddChild(imaTracks);
            numTracks70to90 = numTracks70to90 + 1;
        elseif trackSEL(iTrack,3) >= 0.5*nTimepoints
            imaTracks.SetColor(single(0),single(0),single(1),single(0));
            imaTrackGroup50to70.AddChild(imaTracks);
            numTracks50to70 = numTracks50to70 + 1;
        else
            imaTracks.SetColor(single(1),single(0),single(1),single(0));
            imaTrackGroup0to50.AddChild(imaTracks);
            numTracks0to50 = numTracks0to50 + 1;
        end
        
    end %(for iTrack = 1 : numTracks)

    %give names to groups of tracks
    imaTrackGroup90to100.mName = ['tracks length 90-100% (' ...
        num2str(numTracks90to100) ')'];
    imaTrackGroup70to90.mName = ['tracks length 70-90% (' ...;
        num2str(numTracks70to90) ')'];
    imaTrackGroup50to70.mName = ['tracks length 50-70% (' ...;
        num2str(numTracks50to70) ')'];
    imaTrackGroup0to50.mName = ['tracks length 0-50% (' ...;
        num2str(numTracks0to50) ')'];

    %add track groups to scene
    imaSurpassScene.AddChild(imaTrackGroup90to100);
    imaSurpassScene.AddChild(imaTrackGroup70to90);
    imaSurpassScene.AddChild(imaTrackGroup50to70);
    imaSurpassScene.AddChild(imaTrackGroup0to50);

end %(if select(1))

%% sisters

if select(2)

    %make spots plotted earlier invisible
    imaSpots.mVisible = 0;

    %get tracks from dataStruct
    tracksFinal = dataStruct.tracks;

    %find track start times, end times and lifetimes
    trackSEL = getTrackSEL(tracksFinal);

    %find gaps in tracks
    gapInfo = findTrackGaps(tracksFinal);

    %get list of paired sisters from dataStruct
    pairList = dataStruct.sisterList(1).trackPairs;

    %determine number of pairs
    numPairs = size(pairList,1);

    %create data container for all sisters
    imaTrackAllSisters = imarisApplication.mFactory.CreateDataContainer;
    imaTrackAllSisters.mName = 'sister pairs';

    %plot paired tracks
    for iPair = 1 : numPairs

        %go over both sisters
        for iSister = 1 : 2

            %get track index
            iTrack = pairList(iPair,iSister);

            %create track object
            imaTracks = imarisApplication.mFactory.CreateTrack;
            if iSister == 1
                imaTracks.mName = [num2str(iPair) '_' num2str(iSister) '  (' num2str(pairList(iPair,3)) ')'];
            else
                imaTracks.mName = [num2str(iPair) '_' num2str(iSister)];
            end

            %get spots belonging to this track, where index is per
            %frame
            spotsIndx = [ones(1,trackSEL(iTrack,1)-1) ...
                tracksFinal(iTrack).tracksFeatIndxCG ...
                ones(1,nTimepoints-trackSEL(iTrack,2))]';

            %locate gaps in this track
            gapsInTrack = gapInfo(gapInfo(:,1)==iTrack,:);

            %calculate cumulative index of spots in order to get spot
            %data from the variable "spots"
            spotsIndx = spotsIndx + nSpotSum(1:end-1);

            %get spot coordinates and sizes(some are wrong and will be corrected
            %in the next couple of steps)
            spotsCoord = spots(spotsIndx,1:3);
            spotSize = pixelSize(1)*2*ones(nTimepoints,1);

            %for frames before track starts, assign position as that at
            %the start. Make spot size 0
            spotsCoord(1:trackSEL(iTrack,1)-1,1) = spotsCoord(trackSEL(iTrack,1),1);
            spotsCoord(1:trackSEL(iTrack,1)-1,2) = spotsCoord(trackSEL(iTrack,1),2);
            spotsCoord(1:trackSEL(iTrack,1)-1,3) = spotsCoord(trackSEL(iTrack,1),3);
            spotSize(1:trackSEL(iTrack,1)-1) = 0;

            %for frames after track ends, assign position as that at
            %the end. Make spot size 0
            spotsCoord(trackSEL(iTrack,2)+1:end,1) = spotsCoord(trackSEL(iTrack,2),1);
            spotsCoord(trackSEL(iTrack,2)+1:end,2) = spotsCoord(trackSEL(iTrack,2),2);
            spotsCoord(trackSEL(iTrack,2)+1:end,3) = spotsCoord(trackSEL(iTrack,2),3);
            spotSize(trackSEL(iTrack,2)+1:nTimepoints) = 0;

            %in frames where there is a gap, use coordinate of last
            %frame where object is detected. Make spot size half that of a
            %detected spot
            for iGap = 1 : size(gapsInTrack,1)
                spotsCoord(gapsInTrack(iGap,3):gapsInTrack(iGap,3)+gapsInTrack(iGap,4)-1,1) = spotsCoord(gapsInTrack(iGap,3)-1,1);
                spotsCoord(gapsInTrack(iGap,3):gapsInTrack(iGap,3)+gapsInTrack(iGap,4)-1,2) = spotsCoord(gapsInTrack(iGap,3)-1,2);
                spotsCoord(gapsInTrack(iGap,3):gapsInTrack(iGap,3)+gapsInTrack(iGap,4)-1,3) = spotsCoord(gapsInTrack(iGap,3)-1,3);
                spotSize(gapsInTrack(iGap,3):gapsInTrack(iGap,3)+gapsInTrack(iGap,4)-1) = pixelSize(1);
            end

            %set spot coordinates in imaris object
            imaSpotsTrack = imarisApplication.mFactory.CreateSpots;
            imaSpotsTrack.Set(single(spotsCoord),...
                single(0:nTimepoints-1),single(spotSize));

            %define track spots
            imaTracks.SetSpots(imaSpotsTrack);

            %define track edges
            imaTracks.SetEdges(single([(0:nTimepoints-2)' (1:nTimepoints-1)']));

            %make track invisible
            imaTracks.mVisible = 0;

            %add track to the data container
            imaTrackAllSisters.AddChild(imaTracks);

        end %(for iSister = 1 : 2)
        
    end %(for iPair = 1 : numPairs)

    %make sisters invisible
    imaTrackAllSisters.mVisible = 0;

    %add sisters to scene
    imaSurpassScene.AddChild(imaTrackAllSisters);
    
end %(if select(2))

%% fitted plane

if select(3)
    %disp('Sorry, plotting plane not implemented yet!');
end

% turn warnings back on
warning(warningState)
