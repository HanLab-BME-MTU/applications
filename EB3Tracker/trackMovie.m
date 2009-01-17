function trackMovie(runInfo,indivTrack,timeRange,roiYX,magCoef,showTracks,showDetect)
% TRACKMOVIE makes a movie of all the tracks in a ROI or of an individual
%
%SYNOPSIS trackMovie(runInfo,indivTrack,timeRange,roiYX,magCoef,showTracks,showDetect)
%
%INPUT  runInfo           : structure containing fields .anDir, which gives
%                           the full path to the roi_x directory (under
%                           which there are /feat and /track directories)
%                           and .imDir, which gives the full path to the
%                           folder containing the images for overlay.
%                           if given as [], program will query user for
%                           roi_x directory.
%       indivTrack        : n-vector containing n track numbers if you want
%                           to make n movies of individual tracks. if given
%                           as [], all tracks will be plotted within the
%                           roiYX (if given), within the timeRange
%       timeRange         : row vector of the form [startFrame endFrame]
%                           indicating time range to plot. if not given or 
%                           given as [], tracks from the whole movie will
%                           be displayed. if indivTrack is given, this
%                           parameter defaults to the frames over which the
%                           track appears, regardless of user input.
%       roiYX             : coordinates of a region-of-interest (closed
%                           polygon) in which tracks should be plotted.
%                           if indivTrack is given, this parameter
%                           defaults to the ROI circumscribing the full
%                           track, regardless of user input. if given as
%                           [], user will be asked to select a ROI.
%       magCoef           : factor used to change size of movie on the
%                           screen as it is being produced. if given as [],
%                           the program will make the movie as large as
%                           possible while preserving the aspect ratio.
%       showTracks        : 1 (default) if tracks should be plotted; 0 if not
%       showDetect        : 1 (default) if feature coordinates should be
%                           plotted (color-coded by time) in the movie; 2
%                           if ALL detected features coordinates should be
%                           plotted (not just the ones accepted by the
%                           tracking step) in the movie; 0 if no feature
%                           coordinates should be plotted.
%
%OUTPUT One or more Quicktime movies and the regions of interest used to
%       generate them


close all
homeDir=pwd;

% get runInfo in correct format
if nargin<1 || isempty(runInfo)
    % if not given as input, ask user for ROI directory
    % assume images directory is at same level
    runInfo.anDir=uigetdir(pwd,'Please select ROI directory');
    homeDir=pwd;
    cd(runInfo.anDir);
    cd('..');
    runInfo.imDir=[pwd filesep 'images'];
    cd(homeDir)
else
    % adjust for OS
    if ~isfield(runInfo,'imDir') || ~isfield(runInfo,'anDir')
        error('--trackMovie: first argument should be a structure with fields imDir and anDir');
    else
        [runInfo.anDir] = formatPath(runInfo.anDir);
        [runInfo.imDir] = formatPath(runInfo.imDir);
    end
end

% get output directory from the user
runInfo.movDir=uigetdir(runInfo.anDir,'Please select OUTPUT directory');
cd(runInfo.movDir)

% load movieInfo (detection result)
featDir  = [runInfo.anDir filesep 'feat'];
if ~isdir(featDir)
    error('--trackMovie: feat directory missing')
else
    if exist([featDir filesep 'movieInfo.mat'])
        load([featDir filesep 'movieInfo.mat'])
    else
        error('--trackMovie: movieInfo missing...')
    end
end

% load tracksFinal (tracking result)
trackDir = [runInfo.anDir filesep 'track'];
if ~isdir(trackDir)
    error('--trackMovie: track directory missing')
else
    [listOfFiles]=searchFiles('.mat',[],trackDir,0);
    if ~isempty(listOfFiles)
        load([listOfFiles{1,2} filesep listOfFiles{1,1}])
        if ~exist('tracksFinal','var')
            error('--trackMovie: tracksFinal missing...');
        end
    else
        error('--trackMovie: tracksFinal missing...');
    end
end

% convert tracksFinal to matrix
[trackedFeatureInfo,trackedFeatureIndx] = convStruct2MatNoMS(tracksFinal);
clear trackedFeatureIndx

%get number of tracks and number of time points
[numTracks,numTimePoints] = size(trackedFeatureInfo);
numTimePoints = numTimePoints/8;

% if no individual track given, check/assign timeRange input
if nargin<2 || isempty(indivTrack)
    indivTrack=[];
else
    % only retain tracking info for chosen track(s)
    trackedFeatureInfo=trackedFeatureInfo(indivTrack,:);
end

% check/assign time range
if isempty(indivTrack)
    % check whether a time range for plotting was input
    if nargin<3 || isempty(timeRange)
        timeRange = [1 numTimePoints];
    else
        if timeRange(1) < 1 || timeRange(2) > numTimePoints
            error('--trackMovie: Problem with timeRange');
        end
    end

else
    % assign timeRange for each individual track, disregarding user input
    % begins 3 frames before track appears and ends 3 frames after
    trackSEL = getTrackSEL(trackedFeatureInfo);
    timeRange=zeros(size(trackSEL,1),2);
    for i=1:size(timeRange,1)
        timeRange(i,1)=max(1,trackSEL(i,1)-3);
        timeRange(i,2)=min(numTimePoints,trackSEL(i,2)+3);
    end
end


% extract coordinates - either all tracks if indiv not selected or just the
% individual tracks
tracksX = trackedFeatureInfo(:,1:8:end);
tracksY = trackedFeatureInfo(:,2:8:end);

if nargin<6 || isempty(showTracks)
    showTracks=1;
end

if nargin<7 || isempty(showDetect)
    showDetect=1;
end


[listOfImages] = searchFiles('.tif',[],runInfo.imDir,0);
for iMovie=1:size(timeRange,1)

    startFrame = timeRange(iMovie,1);
    endFrame   = timeRange(iMovie,2);

    % load startFrame image and get size
    fileNameIm = [char(listOfImages(startFrame,2)) filesep char(listOfImages(startFrame,1))];
    img = double(imread(fileNameIm));
    [imL,imW] = size(img);

    if isempty(indivTrack) % all tracks

        % retain all the tracks for this movie
        temptrackedFeatureInfo=trackedFeatureInfo;
        
        if nargin<4 || isempty(roiYX)
            imscaled=(img-min(img(:)))./(max(img(:))-min(img(:)));
            [BW,xi,yi] = roipoly(imscaled);
            roiYX=[yi xi];
        else
            if size(roiYX,2)~=2
                error('--trackMovie: roiYX should be nx2 matrix of coordinates')
            end
        end

        minY=floor(min(roiYX(:,1)));
        maxY=ceil(max(roiYX(:,1)));
        minX=floor(min(roiYX(:,2)));
        maxX=ceil(max(roiYX(:,2)));

        % name the movie with the following:
        % allTracks_startFrame_endFrame_01 (or 02, 03...depending on how
        % many allTracks movies have been produced in the past)

        count=1;
        
        % initialize strings for start and end frame
        s = length(num2str(size(listOfImages,1)));
        strg = sprintf('%%.%dd',s);
        sFrm = sprintf(strg,startFrame);
        eFrm = sprintf(strg,endFrame);
        
         % initialize string for movie number
        s = length(num2str(99)); 
        strg = sprintf('%%.%dd',s);
        movNum = sprintf(strg,count);

        temp=['allTracks_' sFrm '_' eFrm '_' movNum];
        while exist([runInfo.movDir filesep temp '.mov'],'file')>0
            count=count+1;
            movNum = sprintf(strg,count);
            temp=['allTracks_' sFrm '_' eFrm '_' movNum];
        end
        movieName=temp;



    else % specific track

        % just use the current track for plotting
        temptrackedFeatureInfo=trackedFeatureInfo(iMovie,:);
        
        % make roi just around the track coordinates
        minX = floor(max(1,nanmin(tracksX(iMovie,startFrame:endFrame))-10));
        minY = floor(max(1,nanmin(tracksY(iMovie,startFrame:endFrame))-10));
        maxX = ceil(min(imW,nanmax(tracksX(iMovie,startFrame:endFrame))+10));
        maxY = ceil(min(imL,nanmax(tracksY(iMovie,startFrame:endFrame))+10));
        roiYX=[minY minX; maxY minX; maxY maxX; minY maxX; minY minX];

        % name the movie according to track number
        count=1;
        
        % initialize strings for track number
        s = length(num2str(numTracks));
        strg = sprintf('%%.%dd',s);
        trckNum = sprintf(strg,indivTrack(iMovie));
        
        % initialize strings for start and end frame
        s = length(num2str(size(listOfImages,1)));
        strg = sprintf('%%.%dd',s);
        sFrm = sprintf(strg,startFrame);
        eFrm = sprintf(strg,endFrame);
        
         % initialize string for movie number
        s = length(num2str(99)); 
        strg = sprintf('%%.%dd',s);
        movNum = sprintf(strg,count);

        temp=['track_' trckNum '_' sFrm '_' eFrm '_' movNum];
        while exist([runInfo.movDir filesep temp '.mov'],'file')>0
            count=count+1;
            movNum = sprintf(strg,count);
            temp=['track_' trckNum '_' sFrm '_' eFrm '_' movNum];
        end
        movieName=temp;
        
    end

    % calculate the size of the movie from the screen size and roi dimensions
    % if magCoef is given, use it unless bigger than the maximum allowed by the
    % screen size
    scrsz = get(0,'ScreenSize');
    screenW=scrsz(3);
    screenL=scrsz(4);

    if nargin<5 || isempty(magCoef)
        magCoef=inf;
    end

    maxMagCoefW = (0.8*screenW)/(maxX-minX+1);
    maxMagCoefL = (0.8*screenL)/(maxY-minY+1);

    if magCoef > min([maxMagCoefW; maxMagCoefL])
        calcMagCoef = min([magCoef; maxMagCoefW; maxMagCoefL]);
    else
        calcMagCoef = magCoef;
    end

    movieL = (calcMagCoef*(maxY-minY+1));
    movieW = (calcMagCoef*(maxX-minX+1));

    figure('Position',[round(screenW*(1-movieW/screenW)/2) round(screenL*(1-movieL/screenL)/2) movieW movieL])

    eval(['MakeQTMovie start ', movieName '.mov']);

    frmCount2=1;
    colorOverTime = jet(endFrame-startFrame+1);
    for iFrame=startFrame:endFrame

        img=double(imread([char(listOfImages(iFrame,2)) filesep char(listOfImages(iFrame,1))]));
        if showTracks==1
            % [selectedTracks]=plotTracks2D_EB3(trackedFeatureInfo,timeRange,img,ask4sel,plotCurrentOnly,roiYX,movieInfo)
            plotTracks2D_EB3(temptrackedFeatureInfo,[startFrame endFrame],img,0,iFrame,roiYX,[]);
        else
            imagesc(img(minY:maxY,minX:maxX))
            colormap gray
        end

        if showDetect==1 || showDetect==2
            frmCount1=1;
            for j=startFrame:endFrame
                if showDetect==1 % use coordinates from tracks
                    if isempty(indivTrack)
                        xCoord = tracksX(:,j);
                        yCoord = tracksY(:,j);
                    else
                        xCoord = tracksX(iMovie,j);
                        yCoord = tracksY(iMovie,j);
                    end
                elseif showDetect==2 % use coordinates from detection
                    xCoord = vertcat(movieInfo(j).xCoord); xCoord = xCoord(:,1);
                    yCoord = vertcat(movieInfo(j).yCoord); yCoord = yCoord(:,1);
                end
                outOfRangeIdx=find(xCoord<minX | xCoord>maxX | yCoord<minY | yCoord>maxY);
                xCoord(outOfRangeIdx) = [];
                yCoord(outOfRangeIdx) = [];
                xCoord=xCoord-minX+1;
                yCoord=yCoord-minY+1;
                hold on
                scatter(xCoord,yCoord,'Marker','.','cData',repmat(colorOverTime(frmCount1,:),[length(xCoord),1]));
                frmCount1=frmCount1+1;
            end
        end

        text(.25,.25,num2str(iFrame),'Color',colorOverTime(frmCount2,:),'FontWeight','bold','HorizontalAlignment','right','Units','inches')
        MakeQTMovie addaxes
        MakeQTMovie('framerate', 5);
        %MakeQTMovie('quality', .7);

        frmCount2=frmCount2+1;

    end
    MakeQTMovie finish
    save([movieName '_roiYX'],'roiYX')

    close all
end
cd(homeDir)
