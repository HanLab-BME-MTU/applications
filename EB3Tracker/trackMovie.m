function trackMovie(runInfo,movieName,frameRange,showTracks,showDetect,magCoef,roiYX,indivTrack)
% TRACKMOVIE makes a movie of all the tracks in a ROI or of an individual
% track

close all

homeDir=pwd;

if showTracks==1 || showDetect==1
    cd(runInfo.anDir)

    if showTracks==1
        load testTrack % result from Khuloud's tracker
    end
    if showDetect==1
        load(['feat' filesep 'movieInfo']) % detection coordinates
    else
        movieInfo=[];
    end
end

if isempty(movieName)
    movieName='movie';
end
if ~isempty(indivTrack)
    movieName=['movie_' num2str(indivTrack)];
end

% get list of files and image size
[listOfImages] = searchFiles('.tif',[],runInfo.imDir,0);
fileNameIm = [char(listOfImages(1,2)) filesep char(listOfImages(1,1))];
img = double(imread(fileNameIm));
[imL,imW] = size(img);

if isempty(indivTrack)

    % get frame range from image stack or user input
    if isempty(frameRange)
        startFrame = 1;
        endFrame = size(listOfImages,1);
    else
        startFrame = frameRange(1);
        endFrame = frameRange(2);
    end

    if showTracks==1
        % convert structure to matrix, get number of time points
        [trackedFeatureInfo,trackedFeatureIndx] = convStruct2MatNoMS(tracksFinal);
        [numTracks,numTimePoints] = size(trackedFeatureInfo);
        numTimePoints = numTimePoints/8;
        endFrame = min(endFrame,numTimePoints);

    else
        numTimePoints = size(listOfImages,1);
        endFrame = min(endFrame,numTimePoints);
    end

    % select ROI limits for movie
    if isempty(roiYX)
        imscaled=(img-min(img(:)))./(max(img(:))-min(img(:)));
        [BW,xi,yi] = roipoly(imscaled);
        roiYX=[yi xi];
    end
    minY=floor(min(roiYX(:,1)));
    maxY=ceil(max(roiYX(:,1)));
    minX=floor(min(roiYX(:,2)));
    maxX=ceil(max(roiYX(:,2)));

else
    % convert structure to matrix, get number of time points
    [trackedFeatureInfo,trackedFeatureIndx] = convStruct2MatNoMS(tracksFinal);
    [numTracks,numTimePoints] = size(trackedFeatureInfo);
    numTimePoints = numTimePoints/8;

    % only retain tracking info for the single track
    trackedFeatureInfo=trackedFeatureInfo(indivTrack,:);

    % frame range is a few frames before and after track exists
    if isempty(frameRange)
        startFrame=max(1,min(find(~isnan(trackedFeatureInfo(1,1:8:end))))-3);
        endFrame=min(numTimePoints,max(find(~isnan(trackedFeatureInfo(1,1:8:end))))+3);
    else
        startFrame = frameRange(1);
        endFrame = frameRange(2);
    end

    % make roi just around the track coordinates
    tracksX = trackedFeatureInfo(:,1:8:end)';
    tracksY = trackedFeatureInfo(:,2:8:end)';
    minX = floor(max(1,nanmin(tracksX)-10));
    minY = floor(max(1,nanmin(tracksY)-10));
    maxX = ceil(min(imW,nanmax(tracksX)+10));
    maxY = ceil(min(imL,nanmax(tracksY)+10));
    roiYX=[minY minX; maxY minX; maxY maxX; minY maxX; minY minX];

end

% calculate the size of the movie from the screen size and roi dimensions
% if magCoef is given, use it unless bigger than the maximum allowed by the
% screen size
scrsz = get(0,'ScreenSize');
screenW=scrsz(3);
screenL=scrsz(4);

if isempty(magCoef)
    magCoef=inf;
end

maxMagCoefW = (0.8*screenW)/(maxX-minX+1);
maxMagCoefL = (0.8*screenL)/(maxY-minY+1);

if magCoef > min([maxMagCoefW; maxMagCoefL])
    calcMagCoef = min([magCoef; maxMagCoefW; maxMagCoefL]);
else
    calcMagCoef = magCoef;
end

movieL = round(calcMagCoef*(maxY-minY+1));
movieW = round(calcMagCoef*(maxX-minX+1));

figure('Position',[round(screenW*(1-movieW/screenW)/2) round(screenL*(1-movieL/screenL)/2) movieW movieL])


eval(['MakeQTMovie start ', movieName '.avi']);

frmCount2=1;
colorOverTime = jet(endFrame-startFrame+1);
for i=startFrame:endFrame

    img=double(imread([char(listOfImages(i,2)) filesep char(listOfImages(i,1))]));
    if showTracks==1
        % plotTracks2D_EB3(trackedFeatureInfo,timeRange,colorTime,markerType,...
        %    indicateSE,newFigure,img,flipXY,ask4sel,minMaxVel,plotCurrentOnly,roiYX)
        %         [selectedTracks] = plotTracks2D_EB3(trackedFeatureInfo,timeRange,...
        %     newFigure,img,flipXY,ask4sel,plotCurrentOnly,roiYX,movieInfo)
        plotTracks2D_EB3(trackedFeatureInfo,[startFrame endFrame],0,img,0,0,i,roiYX,movieInfo);

        if ~isempty(indivTrack)
            frmCount1=1;
            for j=startFrame:endFrame
                xCoord = tracksX;
                yCoord = tracksY;
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
    else
        imagesc(img(minY:maxY,minX:maxX))
        colormap gray
    end

    if showTracks~=1 && showDetect==1
        frmCount1=1;
        for j=startFrame:endFrame
            xCoord = vertcat(movieInfo(j).xCoord); xCoord = xCoord(:,1);
            yCoord = vertcat(movieInfo(j).yCoord); yCoord = yCoord(:,1);
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
    text(10,10,num2str(i),'Color',colorOverTime(frmCount2,:),'FontWeight','bold','HorizontalAlignment','right')
    MakeQTMovie addaxes
    MakeQTMovie('framerate', 5);
    %MakeQTMovie('quality', .7);

    frmCount2=frmCount2+1;

end
MakeQTMovie finish
save([movieName 'roiYX'],'roiYX')

close all

cd(homeDir)

