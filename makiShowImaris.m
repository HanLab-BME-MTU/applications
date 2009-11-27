function makiShowImaris(dataStruct,select,movieType)
%MAKISHOWIMARIS shows mammalina kinetochore data in Imaris
%
% SYNOPSIS: imarisHandle = makiShowImaris(dataStruct,select)
%
% INPUT dataStruct: (opt) data structure as described in makiMakeDataStruct
%                   if empty, guiLoad
%		select: (opt) vector of analysis results to plot. If not
%               inputed or empty, everything available will be plotted.
%               1st entry: tracks; 2nd entry: sisters, 3rd entry: fitted
%               plane. Enter 1 for result to be plotted, 0 for result not
%               to be plotted.
%       movieType: 1 for deltaVision files, 2 for metamorph stacks.
%                  Optional. Default: 1.
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

if nargin < 3 || isempty(movieType)
    movieType = 1;
end

% turn off property reader warning
warningState = warning;
warning off IMARISIMREAD:NOPROPERTYREADER

% reduce amount of typing
dataProperties = dataStruct.dataProperties;
pixelSize = [dataProperties.PIXELSIZE_XY,dataProperties.PIXELSIZE_XY,...
    dataProperties.PIXELSIZE_Z];

% check for cropping
if isempty(dataProperties.crop)
    dataProperties.crop = zeros(2,3);
end
crop = dataProperties.crop(:,1:3);
isCrop = any(crop,1);
crop(1,~isCrop) = 1;
crop(2,~isCrop) = dataProperties.movieSize(find(~isCrop)); %#ok<FNDSB>

% start imaris
imarisApplication = imarisStartNew;

% read version, which is a string of the form 'Imaris 5.0.0 [May 11 2006]'
% do the easy way and read the 8th and 10th character
%imarisVersion = sscanf(imarisApplication.mVersion,'*s %i.%i.%i %*s %*i');
imarisVersion = [0 0];
vstr = imarisApplication.mVersion;
if ~strcmp(vstr(8),'x')
    imarisVersion(1) = str2num(vstr(8));
    imarisVersion(2) = str2num(vstr(10));
else %KJ: needed when the string is of the form 'Imaris x64 6.3.1 [Jun 16 2009]'
    imarisVersion(1) = str2num(vstr(12));
    imarisVersion(2) = str2num(vstr(14));    
end

interface62 = imarisVersion(1) > 5 && imarisVersion(2) > 1;

% load raw movie into imaris. We could do filtered movie, but
% for this, we would have to load frame by frame and do all the
% image properties stuff
switch movieType
    case 1
        imarisApplication.FileOpen(...
            fullfile(dataStruct.rawMoviePath,dataStruct.rawMovieName),...
            'reader=''DeltaVision''');
    case 2
        imarisApplication.FileOpen(...
            fullfile(dataStruct.rawMoviePath,dataStruct.rawMovieName));
        stackData = metaTiffRead(fullfile(dataStruct.rawMoviePath,dataStruct.rawMovieName),[],[],0);
        numRows = stackData(1).height;
        numCols = stackData(1).width;
end

%         % check image properties: image should begin at -0.5 pix.
%         % It would be nice to be able to set the pixelSize to -0.5.
%         % Instead, we have to read the mins and calculate an offset
%         zeroOffsetX = imarisApplication.mDataSet.mExtendMinX + 0.5*pixelSize(1);
%         zeroOffsetY = imarisApplication.mDataSet.mExtendMinY + 0.5*pixelSize(2);
%         zeroOffsetZ = imarisApplication.mDataSet.mExtendMinZ + 0.5*pixelSize(3);

%get image frame origin offset. I don't think there's any need for the
%subtration done by Jonas, since it is compensated for later - KJ
zeroOffsetX = imarisApplication.mDataSet.mExtendMinX;
zeroOffsetY = imarisApplication.mDataSet.mExtendMinY;
zeroOffsetZ = imarisApplication.mDataSet.mExtendMinZ;
zeroOffset = [zeroOffsetX zeroOffsetY zeroOffsetZ];

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

%% detected spots (initCoord)

%get coordinates
initCoord = dataStruct.initCoord;

if ~isempty(initCoord)
    
    % make spots object: X,Y,Z,T,r
    nTimepoints = dataProperties.movieSize(end);
    nSpots = cat(1,initCoord.nSpots);
    spots = zeros(sum(nSpots),5);
    spots(:,5) = pixelSize(1)*2; % radius in micron
    goodTimes = find(nSpots);
    nSpotSum = [0;cumsum(nSpots)];
    
    for t = goodTimes'
        
        %Jonas's old way. Keep for reference (the -1 here, together with the +0.5
        %above, gives effectively -0.5, as done below) - KJ
        %                 %for DV files, subtract one voxel from the coordinates:
        %                 %The center of the first pixel in Imaris is 0, whereas
        %                 %in our coordinates the center of the first pixel is 1.
        %                 allCoordPix = initCoord(t).allCoordPix(:,[2,1,3]) - 1 + ...
        %                     repmat(crop(1,[2,1,3])-1,nSpots(t),1);
        
        %get positions in pixels, with adjustment for cropping and frame
        %origin (-0.5)
        allCoordPix = initCoord(t).allCoordPix(:,[2,1,3]) - 0.5 + ...
            repmat(crop(1,[2,1,3])-1,nSpots(t),1);
        
        %for STK files, invert imaris y-coordinate
        if movieType==2
            allCoordPix(:,2) = numRows - allCoordPix(:,2);
        end
        
        % calculate positions in microns
        spots(nSpotSum(t)+1:nSpotSum(t+1),1:4) = ...
            [(allCoordPix) .* repmat(pixelSize,nSpots(t),1) + ...
            repmat(zeroOffset,nSpots(t),1),(t-1)*ones(nSpots(t),1)];
        
    end
    
    % create spots object and display
    imaSpotsAll = imarisApplication.mFactory.CreateSpots;
    imaSpotsAll.Set(single(spots(:,1:3)),single(spots(:,4)),single(spots(:,5)));
    imaSpotsAll.mName = ['Spots (avg: ' num2str(size(spots,1)/nTimepoints) ' / frame)'];
    imaSpotsAll.SetColor(single(1),single(1),single(1),single(0));
    imaSurpassScene.AddChild(imaSpotsAll);
    
end


%% tracks

if select(1) && ~isempty(dataStruct.tracks)
    
    if interface62
        % not debugged yet
        
        %make spots plotted earlier invisible
        imaSpotsAll.mVisible = 0;
        
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
        
        %create data container for tracks between 50% and 90% of movie
        imaTrackGroup50to90 = imarisApplication.mFactory.CreateDataContainer;
        
        %create data container for tracks shorter than 50% of movie
        imaTrackGroup0to50 = imarisApplication.mFactory.CreateDataContainer;
        
        %initialize to zero number of tracks in each category
        numTracks90to100 = 0;
        numTracks50to90 = 0;
        numTracks0to50 = 0;
        
        %plot unpaired tracks
        for iTrack = 1 : numTracks
            
            % create 'track' object (keep variable names for simplicity). In
            % 6.2+, use spots object and add edges to it
            imaTracks = imarisApplication.mFactory.CreateSpots;

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
            
            %             % create 'track' object (keep variable names for simplicity). In
            %             % 6.2+, use spots object and add edges to it
            %             imaTracks = imarisApplication.mFactory.CreateSpots;
            
            %define track spots
            imaTracks.Set(single(spotsCoord),...
                single(0:nTimepoints-1),single(spotSize));
            
            %define track edges
            imaTracks.SetTrackEdges(single([(0:nTimepoints-2)' (1:nTimepoints-1)']));
            
            %add track to relevant data container
            if trackSEL(iTrack,3) >= 0.9*nTimepoints
                imaTracks.SetColor(single(1),single(0),single(1),single(0));
                imaTrackGroup90to100.AddChild(imaTracks);
                numTracks90to100 = numTracks90to100 + 1;
            elseif trackSEL(iTrack,3) >= 0.5*nTimepoints
                imaTracks.SetColor(single(0),single(1),single(1),single(0));
                imaTrackGroup50to90.AddChild(imaTracks);
                numTracks50to90 = numTracks50to90 + 1;
            else
                imaTracks.SetColor(single(1),single(1),single(0),single(0));
                imaTrackGroup0to50.AddChild(imaTracks);
                numTracks0to50 = numTracks0to50 + 1;
            end
            
        end %(for iTrack = 1 : numTracks)
        
        %give names to groups of tracks
        imaTrackGroup90to100.mName = ['tracks length 90-100% (' ...
            num2str(numTracks90to100) '-pink)'];
        imaTrackGroup50to90.mName = ['tracks length 50-90% (' ...;
            num2str(numTracks50to90) '-cyan)'];
        imaTrackGroup0to50.mName = ['tracks length 0-50% (' ...;
            num2str(numTracks0to50) '-yellow)'];
        
        %add track groups to scene
        imaSurpassScene.AddChild(imaTrackGroup90to100);
        imaSurpassScene.AddChild(imaTrackGroup50to90);
        imaSurpassScene.AddChild(imaTrackGroup0to50);
        
    else
        % Imaris 5.7.1
        
        %make spots plotted earlier invisible
        imaSpotsAll.mVisible = 0;
        
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
        
        %create data container for tracks between 50% and 90% of movie
        imaTrackGroup50to90 = imarisApplication.mFactory.CreateDataContainer;
        
        %create data container for tracks shorter than 50% of movie
        imaTrackGroup0to50 = imarisApplication.mFactory.CreateDataContainer;
        
        %initialize to zero number of tracks in each category
        numTracks90to100 = 0;
        numTracks50to90 = 0;
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
            elseif trackSEL(iTrack,3) >= 0.5*nTimepoints
                imaTracks.SetColor(single(0),single(0),single(1),single(0));
                imaTrackGroup50to90.AddChild(imaTracks);
                numTracks50to90 = numTracks50to90 + 1;
            else
                imaTracks.SetColor(single(1),single(1),single(0),single(0));
                imaTrackGroup0to50.AddChild(imaTracks);
                numTracks0to50 = numTracks0to50 + 1;
            end
            
        end %(for iTrack = 1 : numTracks)
        
        %give names to groups of tracks
        imaTrackGroup90to100.mName = ['tracks length 90-100% (' ...
            num2str(numTracks90to100) ')'];
        imaTrackGroup50to90.mName = ['tracks length 50-90% (' ...;
            num2str(numTracks50to90) ')'];
        imaTrackGroup0to50.mName = ['tracks length 0-50% (' ...;
            num2str(numTracks0to50) ')'];
        
        %add track groups to scene
        imaSurpassScene.AddChild(imaTrackGroup90to100);
        imaSurpassScene.AddChild(imaTrackGroup50to90);
        imaSurpassScene.AddChild(imaTrackGroup0to50);
        
    end
    
end %(if select(1))


%% sisters

if select(2) && ~isempty(dataStruct.sisterList) && ~isempty(dataStruct.sisterList(1).trackPairs)

    if interface62
        % implementation for Imaris 6.2
        
        %make spots plotted earlier invisible
        imaSpotsAll.mVisible = 0;
        
        %get tracks from dataStruct
        tracksFinal = dataStruct.tracks;
        
        %find track start times, end times and lifetimes
        trackSEL = getTrackSEL(tracksFinal);
        
        %get sister list from dataStruct
        sisterList = dataStruct.sisterList;
        
        %determine number of pairs
        numPairs = length(dataStruct.sisterList);
        
        %create data container for all sisters
        imaTrackAllSisters = imarisApplication.mFactory.CreateDataContainer;
        imaTrackAllSisters.mName = ['sister pairs (' num2str(numPairs) ')'];
        
        %plot paired tracks
        for iPair = 1 : numPairs
            
            %create 'track' object. Make spots; keep old variable names
            imaTracks1 = imarisApplication.mFactory.CreateSpots;
            imaTracks2 = imarisApplication.mFactory.CreateSpots;
            
            %name track objects
            imaTracks1.mName = [num2str(iPair) '_' num2str(1) '  (' num2str(sisterList(1).trackPairs(iPair,3)) ')'];
            imaTracks2.mName = [num2str(iPair) '_' num2str(2)];
            
            %find NaNs in coordinates
            nanCoord = isnan(sisterList(iPair).coords1(:,1));
            
            %determine sister pair start time as the first frame where a zero appears
            startTime = find(nanCoord==0,1,'first');
            
            if ~isempty(startTime) %if sister pair is not completely empty
                
                %determine sister pair end time as the last frame where a zero appears
                endTime = find(nanCoord==0,1,'last');
                
                %locate gaps in sister pair
                missIndx = find(nanCoord==1);
                missIndx = missIndx(missIndx > startTime & missIndx < endTime);
                
                %assign spot sizes (2 pixels in available frames, 1 pixel in
                %missing frames, 0 pixels before and after track)
                spotSize = pixelSize(1)*2*ones(nTimepoints,1);
                spotSize(1:startTime-1) = 0;
                spotSize(endTime+1:end) = 0;
                spotSize(missIndx) = spotSize(missIndx)/2;
                
                %for first sister ...
                
                %get spot indices to obtain positions
                iTrack = sisterList(1).trackPairs(iPair,1);
                shift = trackSEL(iTrack,1) - 1;
                spotsIndx = [ones(1,startTime-1) ...
                    tracksFinal(iTrack).tracksFeatIndxCG(startTime-shift:endTime-shift) ...
                    ones(1,nTimepoints-endTime)]';
                
                %calculate cumulative index of spots in order to get spot
                %data from the variable "spots"
                spotsIndx = spotsIndx + nSpotSum(1:end-1);
                
                %get coordinates (some are wrong and will be corrected
                %in the next couple of steps)
                sisterCoord1 = spots(spotsIndx,1:3);
                
                %for frames before track starts, assign position as that at
                %the start.
                sisterCoord1(1:startTime-1,1) = sisterCoord1(startTime,1);
                sisterCoord1(1:startTime-1,2) = sisterCoord1(startTime,2);
                sisterCoord1(1:startTime-1,3) = sisterCoord1(startTime,3);
                
                %for frames after track ends, assign position as that at
                %the end.
                sisterCoord1(endTime+1:end,1) = sisterCoord1(endTime,1);
                sisterCoord1(endTime+1:end,2) = sisterCoord1(endTime,2);
                sisterCoord1(endTime+1:end,3) = sisterCoord1(endTime,3);
                
                %in frames where there is a gap, use coordinate of last
                %frame where object is detected.
                for iMiss = missIndx'
                    sisterCoord1(iMiss,:) = sisterCoord1(iMiss-1,:);
                end
                
                %set spot coordinates in imaris object
                imaTracks1.Set(single(sisterCoord1),...
                    single(0:nTimepoints-1),single(spotSize));
                
                %define track edges
                imaTracks1.SetTrackEdges(single([(0:nTimepoints-2)' (1:nTimepoints-1)']));
                
                %make track invisible
                imaTracks1.mVisible = 0;
                
                %for second sister ...
                
                %get spot indices to obtain its positions
                iTrack = sisterList(1).trackPairs(iPair,2);
                shift = trackSEL(iTrack,1) - 1;
                spotsIndx = [ones(1,startTime-1) ...
                    tracksFinal(iTrack).tracksFeatIndxCG(startTime-shift:endTime-shift) ...
                    ones(1,nTimepoints-endTime)]';
                
                %calculate cumulative index of spots in order to get spot
                %data from the variable "spots"
                spotsIndx = spotsIndx + nSpotSum(1:end-1);
                
                %get coordinates (some are wrong and will be corrected
                %in the next couple of steps)
                sisterCoord2 = spots(spotsIndx,1:3);
                
                %for frames before track starts, assign position as that at
                %the start.
                sisterCoord2(1:startTime-1,1) = sisterCoord2(startTime,1);
                sisterCoord2(1:startTime-1,2) = sisterCoord2(startTime,2);
                sisterCoord2(1:startTime-1,3) = sisterCoord2(startTime,3);
                
                %for frames after track ends, assign position as that at
                %the end.
                sisterCoord2(endTime+1:end,1) = sisterCoord2(endTime,1);
                sisterCoord2(endTime+1:end,2) = sisterCoord2(endTime,2);
                sisterCoord2(endTime+1:end,3) = sisterCoord2(endTime,3);
                
                %in frames where there is a gap, use coordinate of last
                %frame where object is detected.
                for iMiss = missIndx'
                    sisterCoord2(iMiss,:) = sisterCoord2(iMiss-1,:);
                end
                
                %set spot coordinates in imaris object
                imaTracks2.Set(single(sisterCoord2),...
                    single(0:nTimepoints-1),single(spotSize));
                
                %define track edges
                imaTracks2.SetTrackEdges(single([(0:nTimepoints-2)' (1:nTimepoints-1)']));
                
                %make track invisible
                imaTracks2.mVisible = 0;
                
                %add tracks to the data container
                imaTrackAllSisters.AddChild(imaTracks1);
                imaTrackAllSisters.AddChild(imaTracks2);
                
            end %(if ~isempty(startTime))
            
        end %(for iPair = 1 : numPairs)
        
        %make sisters invisible
        imaTrackAllSisters.mVisible = 0;
        
        %add sisters to scene
        imaSurpassScene.AddChild(imaTrackAllSisters);
        
        
    else
        
        % old version
        
        %make spots plotted earlier invisible
        imaSpotsAll.mVisible = 0;
        
        %get tracks from dataStruct
        tracksFinal = dataStruct.tracks;
        
        %find track start times, end times and lifetimes
        trackSEL = getTrackSEL(tracksFinal);
        
        %get sister list from dataStruct
        sisterList = dataStruct.sisterList;
        
        %determine number of pairs
        numPairs = length(dataStruct.sisterList);
        
        %create data container for all sisters
        imaTrackAllSisters = imarisApplication.mFactory.CreateDataContainer;
        imaTrackAllSisters.mName = ['sister pairs (' num2str(numPairs) ')'];
        
        %plot paired tracks
        for iPair = 1 : numPairs

            %get the sister pair's color
            switch mod(iPair,6)
                case 0
                    color1 = 1;
                    color2 = 0;
                    color3 = 0;
                case 1
                    color1 = 0;
                    color2 = 0;
                    color3 = 1;
                case 2
                    color1 = 1;
                    color2 = 1;
                    color3 = 0;
                case 3
                    color1 = 0;
                    color2 = 1;
                    color3 = 1;
                case 4
                    color1 = 1;
                    color2 = 0;
                    color3 = 1;
                case 5
                    color1 = 1;
                    color2 = 1;
                    color3 = 1;
            end

            %create track object
            imaTracks1 = imarisApplication.mFactory.CreateTrack;
            imaTracks2 = imarisApplication.mFactory.CreateTrack;
            
            %name track objects
            imaTracks1.mName = [num2str(iPair) '_' num2str(1) '  (' num2str(sisterList(1).trackPairs(iPair,3)) ')'];
            imaTracks2.mName = [num2str(iPair) '_' num2str(2)];
            
            %find NaNs in coordinates
            nanCoord = isnan(sisterList(iPair).coords1(:,1));
            
            %determine sister pair start time as the first frame where a zero appears
            startTime = find(nanCoord==0,1,'first');
            
            if ~isempty(startTime) %if sister pair is not completely empty
                
                %determine sister pair end time as the last frame where a zero appears
                endTime = find(nanCoord==0,1,'last');
                
                %locate gaps in sister pair
                missIndx = find(nanCoord==1);
                missIndx = missIndx(missIndx > startTime & missIndx < endTime);
                
                %assign spot sizes (2 pixels in available frames, 1 pixel in
                %missing frames, 0 pixels before and after track)
                spotSize = pixelSize(1)*2*ones(nTimepoints,1);
                spotSize(1:startTime-1) = 0;
                spotSize(endTime+1:end) = 0;
                spotSize(missIndx) = spotSize(missIndx)/2;
                
                %for first sister ...
                
                %get spot indices to obtain positions
                iTrack = sisterList(1).trackPairs(iPair,1);
                shift = trackSEL(iTrack,1) - 1;
                spotsIndx = [ones(1,startTime-1) ...
                    tracksFinal(iTrack).tracksFeatIndxCG(startTime-shift:endTime-shift) ...
                    ones(1,nTimepoints-endTime)]';
                
                %calculate cumulative index of spots in order to get spot
                %data from the variable "spots"
                spotsIndx = spotsIndx + nSpotSum(1:end-1);
                
                %get coordinates (some are wrong and will be corrected
                %in the next couple of steps)
                sisterCoord1 = spots(spotsIndx,1:3);
                
                %for frames before track starts, assign position as that at
                %the start.
                sisterCoord1(1:startTime-1,1) = sisterCoord1(startTime,1);
                sisterCoord1(1:startTime-1,2) = sisterCoord1(startTime,2);
                sisterCoord1(1:startTime-1,3) = sisterCoord1(startTime,3);
                
                %for frames after track ends, assign position as that at
                %the end.
                sisterCoord1(endTime+1:end,1) = sisterCoord1(endTime,1);
                sisterCoord1(endTime+1:end,2) = sisterCoord1(endTime,2);
                sisterCoord1(endTime+1:end,3) = sisterCoord1(endTime,3);
                
                %in frames where there is a gap, use coordinate of last
                %frame where object is detected.
                for iMiss = missIndx'
                    sisterCoord1(iMiss,:) = sisterCoord1(iMiss-1,:);
                end
                
                %set spot coordinates in imaris object
                imaSpotsTrack1 = imarisApplication.mFactory.CreateSpots;
                imaSpotsTrack1.Set(single(sisterCoord1),...
                    single(0:nTimepoints-1),single(spotSize));
                
                %define track spots
                imaTracks1.SetSpots(imaSpotsTrack1);
                
                %define track edges
                imaTracks1.SetEdges(single([(0:nTimepoints-2)' (1:nTimepoints-1)']));
                
                %define track edge color
                imaTracks1.SetColor(single(color1),single(color2),single(color3),single(0));

                %make track invisible
                %                 imaTracks1.mVisible = 0;
                
                %for second sister ...
                
                %get spot indices to obtain its positions
                iTrack = sisterList(1).trackPairs(iPair,2);
                shift = trackSEL(iTrack,1) - 1;
                spotsIndx = [ones(1,startTime-1) ...
                    tracksFinal(iTrack).tracksFeatIndxCG(startTime-shift:endTime-shift) ...
                    ones(1,nTimepoints-endTime)]';
                
                %calculate cumulative index of spots in order to get spot
                %data from the variable "spots"
                spotsIndx = spotsIndx + nSpotSum(1:end-1);
                
                %get coordinates (some are wrong and will be corrected
                %in the next couple of steps)
                sisterCoord2 = spots(spotsIndx,1:3);
                
                %for frames before track starts, assign position as that at
                %the start.
                sisterCoord2(1:startTime-1,1) = sisterCoord2(startTime,1);
                sisterCoord2(1:startTime-1,2) = sisterCoord2(startTime,2);
                sisterCoord2(1:startTime-1,3) = sisterCoord2(startTime,3);
                
                %for frames after track ends, assign position as that at
                %the end.
                sisterCoord2(endTime+1:end,1) = sisterCoord2(endTime,1);
                sisterCoord2(endTime+1:end,2) = sisterCoord2(endTime,2);
                sisterCoord2(endTime+1:end,3) = sisterCoord2(endTime,3);
                
                %in frames where there is a gap, use coordinate of last
                %frame where object is detected.
                for iMiss = missIndx'
                    sisterCoord2(iMiss,:) = sisterCoord2(iMiss-1,:);
                end
                
                %set spot coordinates in imaris object
                imaSpotsTrack2 = imarisApplication.mFactory.CreateSpots;
                imaSpotsTrack2.Set(single(sisterCoord2),...
                    single(0:nTimepoints-1),single(spotSize));
                
                %define track spots
                imaTracks2.SetSpots(imaSpotsTrack2);
                
                %define track edges
                imaTracks2.SetEdges(single([(0:nTimepoints-2)' (1:nTimepoints-1)']));
                
                %define track edge color
                imaTracks2.SetColor(single(color1),single(color2),single(color3),single(0));

                %make track invisible
                %                 imaTracks2.mVisible = 0;
                
                %add tracks to the data container
                imaTrackAllSisters.AddChild(imaTracks1);
                imaTrackAllSisters.AddChild(imaTracks2);
                
            end %(if ~isempty(startTime))
            
        end %(for iPair = 1 : numPairs)
        
        %make sisters invisible
        imaTrackAllSisters.mVisible = 0;
        
        %add sisters to scene
        imaSurpassScene.AddChild(imaTrackAllSisters);
        
    end % switch version
    
end %(if select(2))


%% fitted plane
if select(3) && ~isempty(dataStruct.planeFit)
    
    %make spots plotted earlier invisible
    imaSpotsAll.mVisible = 0;
    
    %get plane fit results and updated classification
    planeFit = dataStruct.planeFit;
    updatedClass = dataStruct.updatedClass;
    
    %generate inlier, unaligned, and lagging indices
    %also get frame phase
    inlierIdx = [];
    unalignedIdx = [];
    laggingIdx = [];
    if isempty(updatedClass)
        for t = goodTimes'
            inlierIdx = [inlierIdx (planeFit(t).inlierIdx + nSpotSum(t))];
            unalignedIdx = [unalignedIdx (planeFit(t).unalignedIdx + nSpotSum(t))];
            laggingIdx = [laggingIdx (planeFit(t).laggingIdx + nSpotSum(t))];
        end
        framePhase = vertcat(planeFit.phase);
    else
        for t = goodTimes'
            inlierIdx = [inlierIdx (updatedClass(t).inlierIdx + nSpotSum(t))];
            unalignedIdx = [unalignedIdx (updatedClass(t).unalignedIdx + nSpotSum(t))];
            laggingIdx = [laggingIdx (updatedClass(t).laggingIdx + nSpotSum(t))];
        end
        framePhase = vertcat(updatedClass.phase);
    end
    
    %find frame where anaphase starts
    aStartFrame = find(framePhase=='a',1,'first');
    
    %plot inliers
    if ~isempty(inlierIdx)
        imaSpotsInlier = imarisApplication.mFactory.CreateSpots;
        imaSpotsInlier.Set(single(spots(inlierIdx,1:3)),single(spots(inlierIdx,4)),single(spots(inlierIdx,5)));
        imaSpotsInlier.mName = ['Inlier spots (avg: ' num2str(length(inlierIdx)/nTimepoints) ' / frame)'];
        imaSpotsInlier.SetColor(single(1),single(1),single(1),single(0));
        imaSurpassScene.AddChild(imaSpotsInlier);
    end
    
    %plot unaligned kinetochores
    if ~isempty(unalignedIdx)
        imaSpotsUnaligned = imarisApplication.mFactory.CreateSpots;
        imaSpotsUnaligned.Set(single(spots(unalignedIdx,1:3)),single(spots(unalignedIdx,4)),...
            single(spots(unalignedIdx,5)));
        imaSpotsUnaligned.mName = ['Unaligned spots (avg: ' num2str(length(unalignedIdx)/nTimepoints) ' / frame)'];
        imaSpotsUnaligned.SetColor(single(1),single(0),single(1),single(0));
        imaSurpassScene.AddChild(imaSpotsUnaligned);
    end
    
    %plot lagging kinetochores
    if ~isempty(laggingIdx)
        imaSpotsLagging = imarisApplication.mFactory.CreateSpots;
        imaSpotsLagging.Set(single(spots(laggingIdx,1:3)),single(spots(laggingIdx,4)),single(spots(laggingIdx,5)));
        imaSpotsLagging.mName = ['Lagging spots (avg: ' num2str(length(laggingIdx)/nTimepoints) ' / frame)'];
        imaSpotsLagging.SetColor(single(0),single(1),single(1),single(0));
        imaSurpassScene.AddChild(imaSpotsLagging);
    end
    
    %initialize variables storing grid spots
    spotsGridEigen = []; %planes from eigenvectors and eigenvalues
    spotsGridInter = []; %planes from interpolation
    
    %go over all time points
    for iTime = 1 : nTimepoints
        
        %if there is a plane to plot ...
        if ~isempty(planeFit(iTime).planeVectors)
            
            %fetch spots in this frame
            spotsFrame = spots(nSpotSum(iTime)+1:nSpotSum(iTime+1),1:3);
            
            %get the origin of the plane to be plotted
            planeOrigin = mean(spotsFrame(planeFit(iTime).inlierIdx,:));
            
            %make the grid
            gridOrigin = planeOrigin([2 1 3]);
            vec1 = 0.2*dataStruct.planeFit(iTime).planeVectors(:,1);
            vec2 = 0.2*dataStruct.planeFit(iTime).planeVectors(:,2);
            vec3 = 0.2*dataStruct.planeFit(iTime).planeVectors(:,3);
            if movieType == 2
                vec1(2) = -vec1(2);
                vec2(2) = -vec2(2);
                vec3(2) = -vec3(2);
            end
            [xGrid,yGrid,zGrid] = arbitraryGrid(vec1,vec2,vec3,...
                gridOrigin,[0 0],[-25 25],[-25 25]);
            
            %assign grid spots coordinates
            spotsGrid = [yGrid(:) xGrid(:) zGrid(:)];
            
            %append frame number and spot size
            spotsGrid = [spotsGrid (iTime-1)*ones(size(spotsGrid,1),1) ...
                pixelSize(1)*0.5*ones(size(spotsGrid,1),1)];
            
            %add grid information to the information on all grids
            if planeFit(iTime).planeVectorClassifier == 1
                spotsGridEigen = [spotsGridEigen; spotsGrid];
            else
                spotsGridInter = [spotsGridInter; spotsGrid];
            end
            
        end %(if ~isempty(planeFit(iTime).planeVectors))
        
    end %(for iTime = 1 : nTimepoints)
    
    %add grids to surpassScene
    if ~isempty(spotsGridEigen)
        imaSpotsGrid = imarisApplication.mFactory.CreateSpots;
        imaSpotsGrid.Set(single(spotsGridEigen(:,1:3)),single(spotsGridEigen(:,4)),single(spotsGridEigen(:,5)));
        imaSpotsGrid.SetColor(single(1),single(1),single(0),single(0));
        if isempty(aStartFrame)
            anaText = [];
        else
            anaText = [' - A: TP ' num2str(aStartFrame)];
        end
        imaSpotsGrid.mName = ['Plane fits (anisotropy)' anaText];
        imaSurpassScene.AddChild(imaSpotsGrid);
    end
    if ~isempty(spotsGridInter)
        imaSpotsGrid2 = imarisApplication.mFactory.CreateSpots;
        imaSpotsGrid2.Set(single(spotsGridInter(:,1:3)),single(spotsGridInter(:,4)),single(spotsGridInter(:,5)));
        imaSpotsGrid2.SetColor(single(1),single(0.75),single(0),single(0));
        imaSpotsGrid2.mName = 'Plane fits (interpolation)';
        imaSurpassScene.AddChild(imaSpotsGrid2);
    end
    
end


% turn warnings back on
warning(warningState)
