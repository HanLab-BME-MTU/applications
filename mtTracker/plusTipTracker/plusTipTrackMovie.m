function [calcMagCoef]=plusTipTrackMovie(projData,indivTrack,timeRange,roiYX,magCoef,showTracks,showDetect,aviInstead,rawToo)
% plusTipTrackMovie makes a movie of all the tracks in a ROI or of an individual
%
%INPUT  projData          : output of plusTipPostTracking, stored in /meta
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
%       rawToo            : 1  to show raw image on left of overlay
%                          {0} to make overlay without dual panel
%
%OUTPUT One or more movies and the regions of interest used to
%       generate them
%       calcMagCoef : magnification factor used to make each movie


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
if ~isfield(projData,'saveDir')
    projData.saveDir=uigetdir(projData.anDir,'Please select OUTPUT directory');
end
if isequal(projData.saveDir,0)
    disp('No output directory selected.')
    return
end
cd(projData.saveDir)


trackData=projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix;
% if no individual track given, use all

% assign frame range
if nargin<3 || isempty(timeRange)
    timeRange=[1 projData.nFrames];
else
    if timeRange(1)<1
        timeRange(1)=1;
    end
    if timeRange(2)>projData.nFrames
        timeRange(2)=projData.nFrames;
    end
end

if nargin<2 || isempty(indivTrack)
    indivTrack=[];
    subIdx{1,1}=[];

    % get coordinates of track
    tracksX = projData.xCoord;
    tracksY = projData.yCoord;


else % otherwise, use individul tracks
    timeRangeIndiv=zeros(length(indivTrack),2);
    for i=1:length(indivTrack)
        % cell array with subtrack indices for each track
        subIdx{i,1}=find(trackData(:,1)==indivTrack(i));
        % start and end frames for whole track
        sF=min(trackData(subIdx{i,1},2));
        eF=max(trackData(subIdx{i,1},3));
        % assign frame range for each individual track, bounded by user
        % input
        % begins 3 frames before track appears and ends 3 frames after
        timeRangeIndiv(i,1)=max(timeRange(1),sF-3); % could do sF-3 and eF+3
        timeRangeIndiv(i,2)=min(timeRange(2),eF+3);
        if  timeRangeIndiv(i,1)>=timeRangeIndiv(i,2)
            error(['--plusTipTrackMovie: check range for individual track ' num2str(indivTrack(i))]);
        end
    end
    clear timeRange
    timeRange=timeRangeIndiv;

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

if aviInstead==1
    movStr='.avi';
else
    movStr='.mov';
end

% load startFrame image and get size
[listOfImages] = searchFiles('.tif',[],projData.imDir,0);
nImTot = size(listOfImages,1);

imageName1 = [char(listOfImages(1,2)) filesep char(listOfImages(1,1))];
 

%Check if Image Name Numbers are Padded (Ext no will be greater than
%1 if padded) If Image Name Not Padded Than Sort

[path body no ext ] = getFilenameBody(imageName1);

if length(no)>1
   fileNameIm = [char(listOfImages(1,2)) filesep char(listOfImages(1,1))];
   padded = 1;
else 
    padded = 0;
end

if padded == 0     
 %Initialize Cells for Sorting
    %path = path of file, body = body of filename, ext = extension of filename 
    %(tif etc) (all of these require a cell because they are strings)
    % num = number of filename (do not want in cell so can sort)
    pathCell = cell(nImTot,1);
    bodyCell = cell(nImTot,1);
    extCell = cell(nImTot,1);
    num = zeros(nImTot,1);

        %Sort List
        % For each frame get the image name from listOfImages
        for iFrame =  1:nImTot;
            imageName = [char(listOfImages(iFrame,2)) filesep char(listOfImages(iFrame,1))];

            %Call "getFilenameBody" from common dir to split filename into path,
            %body, no, and ext. Put path, body and ext into appropriate cell/number vector for that
            %frame
            [path body no ext ] = getFilenameBody(imageName);


            pathCell(iFrame) = cellstr(path);
            bodyCell(iFrame) = cellstr(body);
            extCell(iFrame) = cellstr(ext);

            % output "no" is a string so convert to number to sort
            num(iFrame)  = str2double(no);
 
        end
    %Sort number vector numerically
    sortednum = sort(num);

    %Convert number vector to cell
    sortednum_cell = num2cell(sortednum);

    %Create Sorted Image List
    sortedImages = [pathCell, bodyCell, sortednum_cell, extCell];
    
    fileNameIm = [char(sortedImages(1,1)) filesep char(sortedImages(1,2)),...;
    num2str(sortednum(1)) char(sortedImages(1,4))];
end


img1 = double(imread(fileNameIm));
[imL,imW,dim] = size(img1);

calcMagCoef=magCoef;
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
        elseif islogical(roiYX)
            [r c]=find(roiYX);
            roiYX=[min(r) min(c); max(r) max(c)];
        else
            if size(roiYX,2)~=2
                error('--plusTipTrackMovie: roiYX should be nx2 matrix of coordinates or logical mask')
            end
        end


        minY=floor(min(roiYX(:,1)));
        maxY=ceil(max(roiYX(:,1)));
        minX=floor(min(roiYX(:,2)));
        maxX=ceil(max(roiYX(:,2)));

        % name the movie with the following:
        % allTracks_startFrame_endFrame_01 (or 02, 03...depending on how
        % many allTracks movies have been produced in the past)
        temp=['allTracks_' sFrm '_' eFrm '_' movNum];
        while exist([projData.saveDir filesep temp movStr],'file')>0
            count=count+1;
            movNum = sprintf(strg,count);
            temp=['allTracks_' sFrm '_' eFrm '_' movNum];
        end
        movieName=temp;

    else % specific track

        [imL,imW,dim] = size(img1);
        % make roi just around the track coordinates with 10pix cushion
        minX = floor(max(1,nanmin(tracksX(iMovie,startFrame:endFrame))-10));
        minY = floor(max(1,nanmin(tracksY(iMovie,startFrame:endFrame))-10));
        maxX = ceil(min(imW,nanmax(tracksX(iMovie,startFrame:endFrame))+10));
        maxY = ceil(min(imL,nanmax(tracksY(iMovie,startFrame:endFrame))+10));
        roiYX=[minY minX; maxY maxX];

        % initialize strings for track number
        s = length(num2str(projData.nTracks));
        strg1 = sprintf('%%.%dd',s);
        trckNum = sprintf(strg1,indivTrack(iMovie));

        temp=['track_' trckNum '_' sFrm '_' eFrm '_' movNum];
        while exist([projData.saveDir filesep temp movStr],'file')>0
            count=count+1;
            movNum = sprintf(strg,count);
            temp=['track_' trckNum '_' sFrm '_' eFrm '_' movNum];
        end
        movieName=temp;

    end

    % calculate the size of the movie from the screen size and roi dimensions
    % if magCoef is given, use it unless bigger than the maximum allowed by the
    % screen size
    if rawToo==1
        xRange=2*(maxX-minX+1);
        yRange=maxY-minY+1;
    else
        xRange=maxX-minX+1;
        yRange=maxY-minY+1;
    end

    maxMagCoefW = (0.8*screenW)/xRange;
    maxMagCoefL = (0.8*screenL)/yRange;

    if magCoef(iMovie) > min([maxMagCoefW; maxMagCoefL])
        calcMagCoef(iMovie) = min([magCoef(iMovie); maxMagCoefW; maxMagCoefL]);
    else
        calcMagCoef(iMovie) = magCoef(iMovie);
    end

    movieL = (calcMagCoef(iMovie)*yRange);
    movieW = (calcMagCoef(iMovie)*xRange);

    figure('Position',[round(screenW*(1-movieW/screenW)/2) round(screenL*(1-movieL/screenL)/2) movieW movieL])

    % if QT, initiate start of the movie
    if aviInstead==0
        eval(['MakeQTMovie start ', movieName '.mov']);
    end

    frmCount=1;
    colorOverTime = jet(endFrame-startFrame+1);
    for iFrame=startFrame:endFrame
        if padded == 1 % if file number padded with zeros use original file list
            img=(imread([char(listOfImages(iFrame,2)) filesep char(listOfImages(iFrame,1))]));
        else % use sorted list
            img = (imread([char(sortedImages(iFrame,1)) filesep char(sortedImages(iFrame,2)),...;
            num2str(sortednum(iFrame)) char(sortedImages(iFrame,4))]));
        end 
        if showTracks==1
            % plot the tracks
            plusTipPlotTracks(projData,subIdx{iMovie,1},[startFrame endFrame],img,0,iFrame,roiYX,[],rawToo);
            % crop the image based on roi but avoid making data type double
            % also make sure it works for BW or RGB images
            [imL,imW,dim]=size(img);
            temp=zeros(imL,imW);
            temp(minY:maxY,minX:maxX)=1;
            idx=find(temp(:));
            idxAll=[];
            for i=1:dim
                idxAll=[idxAll; idx+(i-1)*(imL*imW)];
            end
            img=reshape(img(idxAll),maxY-minY+1,maxX-minX+1,[]);

            if rawToo==1
                img=[img img];
            end

        else
            % otherwise just show the image
            % crop the image based on roi but avoid making data type double
            % also make sure it works for BW or RGB images
            [imL,imW,dim]=size(img);
            temp=zeros(imL,imW);
            temp(minY:maxY,minX:maxX)=1;
            idx=find(temp(:));
            idxAll=[];
            for i=1:dim
                idxAll=[idxAll; idx+(i-1)*(imL*imW)];
            end
            img=reshape(img(idxAll),maxY-minY+1,maxX-minX+1,[]);
            if rawToo==1
                img=[img img];
            end
            if size(img,3)==1
                imagesc(double(img));
                colormap gray
            else
                image(img);
            end
            %axis equal
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
                    try
                        xCoord=arrayfun(@(i) vertcat(movieInfo(i).xCoord(:,1)),[startFrame:endFrame]','UniformOutput',0);
                        yCoord=arrayfun(@(i) vertcat(movieInfo(i).yCoord(:,1)),[startFrame:endFrame]','UniformOutput',0);
                        l=arrayfun(@(i) length(xCoord{i}),[1:nFr]');
                        col=cell2mat(arrayfun(@(i) repmat(colorOverTime(i,:),[l(i) 1]), [1:nFr]','UniformOutput',0));
                        xCoord=cell2mat(xCoord);
                        yCoord=cell2mat(yCoord);
                    catch
                        xCoord=[];
                        yCoord=[];
                        l=[];
                        col=[];
                    end

                case 3 % plot all features detected in current frame only in dark blue
                    nFr=length(iFrame:iFrame);
                    try % in case no features in current frame
                        xCoord=arrayfun(@(i) vertcat(movieInfo(i).xCoord(:,1)),[iFrame:iFrame]','UniformOutput',0);
                        yCoord=arrayfun(@(i) vertcat(movieInfo(i).yCoord(:,1)),[iFrame:iFrame]','UniformOutput',0);
                        l=arrayfun(@(i) length(xCoord{i}),[1:nFr]');
                        col=cell2mat(arrayfun(@(i) repmat(colorOverTime(1,:),[l(i) 1]), [1:nFr]','UniformOutput',0));
                        xCoord=cell2mat(xCoord);
                        yCoord=cell2mat(yCoord);
                    catch
                        xCoord=[];
                        yCoord=[];
                        l=[];
                        col=[];
                    end

            end

            % remove those coordinates out of the ROI
            outOfRangeIdx=find(xCoord<minX | xCoord>maxX | yCoord<minY | yCoord>maxY);
            xCoord(outOfRangeIdx) = [];
            yCoord(outOfRangeIdx) = [];
            col(outOfRangeIdx,:) = [];
            xCoord=xCoord-minX+1;
            yCoord=yCoord-minY+1;

            if rawToo==1
                xCoord=xCoord+length(img(1,:,1))/2;
                % y values don't change
            end

            % pot detection points
            hold on
            scatter(xCoord,yCoord,'Marker','.','cData',col);

        end

        % plot color-coded frame number
        text(.25,.25,num2str(iFrame),'Color',colorOverTime(frmCount,:),'FontWeight','bold','HorizontalAlignment','right','Units','inches')

        % make dark line down the middle of split-panel
        if rawToo==1
            [imL,imW,dim]=size(img);
            plot([imW/2; imW/2],[0; imL+1],'k','lineWidth',3)
        end

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
        %movie2avi       (F,[projData.saveDir filesep movieName '.avi'],'COMPRESSION','Cinepak','FPS',5)
        movie2aviNADA_CAW(F,[projData.saveDir filesep movieName '.avi'],'COMPRESSION','Cinepak','FPS',5)
    else
        MakeQTMovie finish
    end

    save([movieName '_roiYX'],'roiYX')

    close(figure(gcf));
    clear F
end
cd(homeDir)

