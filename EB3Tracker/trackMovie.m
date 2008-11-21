function trackMovie(runInfo,movieName,frameRange,showTracks,showDetect,magCoef,roiYX)
% TRACKMOVIE makes a movie of all the tracks in a ROI

close all

homeDir=pwd;

if isempty(runInfo)
    runInfo.imDir = uigetdir(pwd,'Select image directory');
    if showTracks==1 || showDetect==1
        cd(runInfo.imDir)
        cd ..
        cd 'roi_analysis'
        cd feat
    end
    if showTracks==1
        load testTrack050
    end
    if showDetect==1
        load movieInfo
    else
        movieInfo=[];
    end
else
    if showTracks==1 || showDetect==1
        cd([runInfo.anDir filesep 'feat'])
    end
    if showTracks==1
        load testTrack050_2
    end
    if showDetect==1
        load movieInfo
    else
        movieInfo=[];
    end
end

if isempty(movieName)
    movieName='movie';
end


% get list of files and image size
[listOfImages] = searchFiles('.tif',[],runInfo.imDir,0);
fileNameIm = [char(listOfImages(1,2)) filesep char(listOfImages(1,1))];
img = double(imread(fileNameIm));
[imL,imW] = size(img);

if isempty(frameRange)
    startFrame = 1;
    endFrame = size(listOfImages,1);
else
    startFrame = frameRange(1);
    endFrame = frameRange(2);
end

% select ROI limits for movie
if nargin<7
imscaled=(img-min(img(:)))./(max(img(:))-min(img(:)));
[BW,xi,yi] = roipoly(imscaled);
roiYX=[yi xi];
end
minY=floor(min(roiYX(:,1)));
maxY=ceil(max(roiYX(:,1)));
minX=floor(min(roiYX(:,2)));
maxX=ceil(max(roiYX(:,2)));

close all

scrsz = get(0,'ScreenSize');
screenW=scrsz(3);
screenL=scrsz(4);
movieL = magCoef*(maxY-minY+1);
movieW = magCoef*(maxX-minX+1);

figure('Position',[round(screenW*(1-movieW/screenW)/2) round(screenL*(1-movieL/screenL)/2) movieW movieL])

if showTracks==1
    % convert tracksFinal structure to matrix
    [trackedFeatureInfo,trackedFeatureIndx] = convStruct2MatNoMS(tracksFinal);
    [numTracks,numTimePoints] = size(trackedFeatureInfo);
    numTimePoints = numTimePoints/8;
    endFrame = min(endFrame,numTimePoints);

else
    numTimePoints = size(listOfImages,1);
    endFrame = min(endFrame,numTimePoints);
end

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
    text(15,15,num2str(i),'Color',colorOverTime(frmCount2,:),'FontWeight','bold','HorizontalAlignment','right')
    MakeQTMovie addaxes
    MakeQTMovie('framerate', 5);
    %MakeQTMovie('quality', .7);
    
    frmCount2=frmCount2+1;

end
MakeQTMovie finish

close all

cd(homeDir)

