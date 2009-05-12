function plusTipTrackMovie(projData,indivTrack,timeRange,roiYX,magCoef,showTracks,showDetect,aviInstead)
% plusTipTrackMovie makes a movie of all the tracks in a ROI or of an individual
%
%INPUT  projData          : output of metaEB3analysis, stored in /meta
%                           if not given, user will be asked to select the
%                           file
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
%                           polygon) in which tracks should be plotted, or
%                           a logical mask the size of the image. if
%                           indivTrack is given, this parameter
%                           defaults to the ROI circumscribing the full
%                           track, regardless of user input. if given as
%                           [], user will be asked to select a ROI.
%       magCoef           : factor used to change size of movie on the
%                           screen as it is being produced. if given as [],
%                           the program will make the movie as large as
%                           possible while preserving the aspect ratio.
%       showTracks        : 0  tracks not plotted
%                          {1} tracks plotted
%       showDetect        : 0  no features centroids plotted
%                          {1} ONLY feature centroids used in tracks plotted,
%                              color-coded by frame
%                           2  ALL detected feature centroids plotted,
%                              color-coded by frame
%                           3  ALL detected feature centroids plotted,
%                              current frame only (dark blue)
%       aviInstead        : 1  to make AVI movie (works in Windows only),
%                          {0} to make MOV movie
%
%OUTPUT One or more movies and the regions of interest used to
%       generate them


homeDir=pwd;

% get projData in correct format
if nargin<1 || isempty(projData)
    % if not given as input, ask user for ROI directory
    % assume images directory is at same level
    [fileName,pathName]=uigetfile('*.mat','Please select projData from META directory');
    if isequal(fileName,0)
        return
    end
    projData=load([pathName filesep fileName]);
    projData=projData.projData;
end
projData.anDir=formatPath(projData.anDir);
projData.imDir=formatPath(projData.imDir);

% get output directory from the user
projData.movDir=uigetdir(projData.anDir,'Please select OUTPUT directory');
if isequal(projData.movDir,0)
    return
end
cd(projData.movDir)


trackData=projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix;
% if no individual track given, use all
if nargin<2 || isempty(indivTrack)
    indivTrack=[];
    subIdx{1,1}=[];

    % assign frame range
    if nargin<3 || isempty(timeRange)
        timeRange=[1 projData.numFrames];
    else
        if timeRange(1)<1
            timeRange(1)=1;
        end
        if timeRange(2)>projData.numFrames
            timeRange(2)=projData.numFrames;
        end
    end

    % get coordinates of track
    tracksX = projData.xCoord;
    tracksY = projData.yCoord;


else % otherwise, use individul tracks
    timeRange=zeros(length(indivTrack),2);
    for i=1:length(indivTrack)
        % cell array with subtrack indices for each track
        subIdx{i,1}=find(trackData(:,1)==indivTrack(i));
        % start and end frames for whole track
        sF=min(trackData(subIdx{i,1},2));
        eF=max(trackData(subIdx{i,1},3));
        % assign frame range for each individual track, disregarding user input
        % begins 3 frames before track appears and ends 3 frames after
        timeRange(i,1)=max(1,sF-3);
        timeRange(i,2)=min(projData.numFrames,eF+3);
    end

    % get coordinates of track
    tracksX = projData.xCoord(indivTrack,:);
    tracksY = projData.yCoord(indivTrack,:);

end

% check input for magnification coefficient
if nargin<5 || isempty(magCoef)
    magCoef=inf;
end
scrsz = get(0,'ScreenSize');
screenW=scrsz(3);
screenL=scrsz(4);


if nargin<6 || isempty(showTracks)
    showTracks=1;
end

if nargin<7 || isempty(showDetect)
    showDetect=1;
end
if showDetect==2 || showDetect==3
    % load movieInfo (detection result)
    featDir  = [projData.anDir filesep 'feat'];
    if ~isdir(featDir)
        error('--plusTipTrackMovie: feat directory missing')
    else
        if exist([featDir filesep 'movieInfo.mat'])
            load([featDir filesep 'movieInfo.mat'])
        else
            error('--plusTipTrackMovie: movieInfo missing...')
        end
    end
end

if nargin<8 || isempty(aviInstead)
    aviInstead=0;
end

% load startFrame image and get size
[listOfImages] = searchFiles('.tif',[],projData.imDir,0);
fileNameIm = [char(listOfImages(timeRange(1,1),2)) filesep char(listOfImages(timeRange(1,1),1))];
img1 = double(imread(fileNameIm));
[imL,imW] = size(img1);

for iMovie=1:size(timeRange,1)

    startFrame = timeRange(iMovie,1);
    endFrame   = timeRange(iMovie,2);

    % count used to make multiple movies with same prefix in name
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



    if isempty(indivTrack) % all tracks
        % check input for ROI coordinates
        if nargin<4 || isempty(roiYX)
            imscaled=(img1-min(img1(:)))./(max(img1(:))-min(img1(:)));
            figure
            set(gcf,'Name','Choose ROI by clicking on image...')
            [BW,xi,yi] = roipoly(imscaled);
            close(figure(gcf))
            roiYX=[yi xi];
        end

        minY=floor(min(roiYX(:,1)));
        maxY=ceil(max(roiYX(:,1)));
        minX=floor(min(roiYX(:,2)));
        maxX=ceil(max(roiYX(:,2)));

        % name the movie with the following:
        % allTracks_startFrame_endFrame_01 (or 02, 03...depending on how
        % many allTracks movies have been produced in the past)
        temp=['allTracks_' sFrm '_' eFrm '_' movNum];
        while exist([projData.movDir filesep temp '.mov'],'file')>0
            count=count+1;
            movNum = sprintf(strg,count);
            temp=['allTracks_' sFrm '_' eFrm '_' movNum];
        end
        movieName=temp;

    else % specific track

        % make roi just around the track coordinates with 10pix cushion
        minX = floor(max(1,nanmin(tracksX(iMovie,startFrame:endFrame))-10));
        minY = floor(max(1,nanmin(tracksY(iMovie,startFrame:endFrame))-10));
        maxX = ceil(min(imW,nanmax(tracksX(iMovie,startFrame:endFrame))+10));
        maxY = ceil(min(imL,nanmax(tracksY(iMovie,startFrame:endFrame))+10));
        roiYX=[minY minX; maxY maxX];

        % initialize strings for track number
        s = length(num2str(projData.numTracks));
        strg = sprintf('%%.%dd',s);
        trckNum = sprintf(strg,indivTrack(iMovie));

        temp=['track_' trckNum '_' sFrm '_' eFrm '_' movNum];
        while exist([projData.movDir filesep temp '.mov'],'file')>0
            count=count+1;
            movNum = sprintf(strg,count);
            temp=['track_' trckNum '_' sFrm '_' eFrm '_' movNum];
        end
        movieName=temp;

    end

    % calculate the size of the movie from the screen size and roi dimensions
    % if magCoef is given, use it unless bigger than the maximum allowed by the
    % screen size
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

    % if QT, initiate start of the movie
    if aviInstead==0
        eval(['MakeQTMovie start ', movieName '.mov']);
    end

    frmCount=1;
    colorOverTime = jet(endFrame-startFrame+1);
    for iFrame=startFrame:endFrame

        img=double(imread([char(listOfImages(iFrame,2)) filesep char(listOfImages(iFrame,1))]));
        if showTracks==1
            % plot the tracks
            plusTipPlotTracks(projData,subIdx{iMovie,1},[startFrame endFrame],img,0,iFrame,roiYX,[]);
        else
            % otherwise just show the image
            imagesc(img(minY:maxY,minX:maxX))
            colormap gray
        end
        if showDetect~=0
            switch showDetect
                case 1 % plot detection features used in tracks, color-coded by frame
                    if isempty(indivTrack)
                        xCoord=tracksX(:,startFrame:endFrame);
                        yCoord=tracksY(:,startFrame:endFrame);
                    else
                        xCoord=tracksX(iMovie,startFrame:endFrame);
                        yCoord=tracksY(iMovie,startFrame:endFrame);
                    end
                    [nTracks,nFr]=size(xCoord);
                    col=cell2mat(arrayfun(@(i) repmat(colorOverTime(i,:),[nTracks 1]), [1:nFr]','UniformOutput',0));
                    xCoord=xCoord(:);
                    yCoord=yCoord(:);

                case 2 % plot all features detected, color-coded by frame
                    nFr=length(startFrame:endFrame);
                    xCoord=arrayfun(@(i) vertcat(movieInfo(i).xCoord(:,1)),[startFrame:endFrame]','UniformOutput',0);
                    yCoord=arrayfun(@(i) vertcat(movieInfo(i).yCoord(:,1)),[startFrame:endFrame]','UniformOutput',0);
                    l=arrayfun(@(i) length(xCoord{i}),[1:nFr]');
                    col=cell2mat(arrayfun(@(i) repmat(colorOverTime(i,:),[l(i) 1]), [1:nFr]','UniformOutput',0));
                    xCoord=cell2mat(xCoord);
                    yCoord=cell2mat(yCoord);

                case 3 % plot all features detected in current frame only in dark blue
                    nFr=length(iFrame:iFrame);
                    xCoord=arrayfun(@(i) vertcat(movieInfo(i).xCoord(:,1)),[iFrame:iFrame]','UniformOutput',0);
                    yCoord=arrayfun(@(i) vertcat(movieInfo(i).yCoord(:,1)),[iFrame:iFrame]','UniformOutput',0);
                    l=arrayfun(@(i) length(xCoord{i}),[1:nFr]');
                    col=cell2mat(arrayfun(@(i) repmat(colorOverTime(1,:),[l(i) 1]), [1:nFr]','UniformOutput',0));
                    xCoord=cell2mat(xCoord);
                    yCoord=cell2mat(yCoord);
            end

            % remove those coordinates out of the ROI
            outOfRangeIdx=find(xCoord<minX | xCoord>maxX | yCoord<minY | yCoord>maxY);
            xCoord(outOfRangeIdx) = [];
            yCoord(outOfRangeIdx) = [];
            col(outOfRangeIdx,:) = [];
            xCoord=xCoord-minX+1;
            yCoord=yCoord-minY+1;

            % pot detection points
            hold on
            scatter(xCoord,yCoord,'Marker','.','cData',col);

        end

        % plot color-coded frame number
        text(.25,.25,num2str(iFrame),'Color',colorOverTime(frmCount,:),'FontWeight','bold','HorizontalAlignment','right','Units','inches')

        % add frame to the movie
        if aviInstead==1
            F(frmCount) = getframe;
        else
            MakeQTMovie addaxes
            MakeQTMovie('framerate', 5);
            %MakeQTMovie('quality', .7);
        end

        frmCount=frmCount+1;

    end

    if aviInstead==1
        %movie2avi       (F,[projData.movDir filesep movieName '.avi'],'COMPRESSION','Cinepak','FPS',5)
        movie2aviNADA_CAW(F,[projData.movDir filesep movieName '.avi'],'COMPRESSION','Cinepak','FPS',5)
    else
        MakeQTMovie finish
    end

    save([movieName '_roiYX'],'roiYX')

    close(figure(gcf));
    clear F
end
cd(homeDir)

