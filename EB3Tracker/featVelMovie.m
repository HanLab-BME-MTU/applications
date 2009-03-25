function featVelMovie(projData,timeRange,velLimit,roiYX,aviInstead)
% FEATVELMOVIE makes a movie of features color-coded by velocity
%
%SYNOPSIS featVelMovie(projData,timeRange,velLimit,roiYX)
%
%INPUT  projData          : output of metaEB3analysis, stored in /meta
%                           if not given, user will be asked to select the
%                           file
%       timeRange         : row vector of the form [startFrame endFrame]
%                           indicating time range to plot. if not given or 
%                           given as [], tracks from the whole movie will
%                           be displayed.
%       velLimit          : max speed to use for colormap (speeds higher
%                           than velLimit or lower than -velLimit will be
%                           plotted as red or blue, respectively.
%       roiYX             : coordinates of a region-of-interest (closed
%                           polygon) in which tracks should be plotted, or
%                           a logical mask the size of the image. if given
%                           as [], user will be asked to select a ROI.
%       aviInstead        : 1 to make AVI move (works in Windows only),
%                           0 (default) for Quicktime
%
%OUTPUT Quicktime movies and the regions of interest used to generate them.
%       Cool colors correspond to shrinkage events (these have "negative"
%       speeds); warm colors correspond to growth or pause.

close all
homeDir=pwd;

% get runInfo in correct format
if nargin<1 || isempty(projData)
    % if not given as input, ask user for ROI directory
    % assume images directory is at same level
    [fileName,pathName]=uigetfile('*.mat','Please select projData from META directory');
    projData=load([pathName filesep fileName]);
    projData=projData.projData;
end


projData.anDir=formatPath(projData.anDir);
projData.imDir=formatPath(projData.imDir);

% get output directory from the user
movDir=[projData.anDir filesep 'movies'];
if ~isdir(movDir)
    mkdir(movDir)
end
cd(movDir)

% check whether a time range for plotting was input
if nargin<3 || isempty(timeRange)
    timeRange = [1 projData.numFrames];
else
    if timeRange(1) < 1 || timeRange(2) > projData.numFrames
        error('--propVelScatter: Problem with timeRange');
    end
end
startFrame = timeRange(1);
endFrame   = timeRange(2);
    

% default: use maximum velocity as limit
if nargin<3 || isempty(velLimit)
    velLimit = inf;
end

% get list of files and image size
[listOfImages] = searchFiles('.tif',[],projData.imDir,0);
fileNameIm = [char(listOfImages(1,2)) filesep char(listOfImages(1,1))];
img = double(imread(fileNameIm));
[imL,imW] = size(img);


% check input for ROI coordinates
if nargin<4 || isempty(roiYX)
    imscaled=(img-min(img(:)))./(max(img(:))-min(img(:)));
    [BW,xi,yi] = roipoly(imscaled);
    roiYX=[yi xi];
elseif islogical(roiYX)
    [r c]=find(roiYX);
    roiYX=[min(r) min(c); max(r) max(c)];
else
    if size(roiYX,2)~=2
        error('--featVelMovie: roiYX should be nx2 matrix of coordinates or logical mask')
    end
end

if nargin<5 || isempty(aviInstead)
    aviInstead=0;
end

minY=floor(min(roiYX(:,1)));
maxY=ceil(max(roiYX(:,1)));
minX=floor(min(roiYX(:,2)));
maxX=ceil(max(roiYX(:,2)));


% get velocity limit for color map
cMapLength=128; cMap=jet(cMapLength);
trackVelMatrix=projData.segGapAvgVel_micPerMin;
m=max(abs([min(trackVelMatrix(:)); max(trackVelMatrix(:))]));
m=round(min(m,velLimit));
trackVelMatrix(trackVelMatrix<-m)=-m;
trackVelMatrix(trackVelMatrix>m)=m;
mapper=linspace(-m,m,cMapLength)';


% name the movie according to track number
count=1;

% initialize strings for start and end frame
s = length(num2str(size(listOfImages,1)));
strg = sprintf('%%.%dd',s);
sFrm = sprintf(strg,startFrame);
eFrm = sprintf(strg,endFrame);

% initialize string for velocity limit
s = length(num2str(99));
strg = sprintf('%%.%dd',s);
velLim = sprintf(strg,m);

% initialize string for movie number
s = length(num2str(99));
strg = sprintf('%%.%dd',s);
movNum = sprintf(strg,count);

% name the movie
temp=['trackVelocityMovie_' sFrm '_' eFrm '_' velLim '_' movNum];
while exist([movDir filesep temp '.mov'],'file')>0
    count=count+1;
    movNum = sprintf(strg,count);
    temp=['trackVelocityMovie_' sFrm '_' eFrm '_' velLim '_' movNum];
end
movieName=temp;


% calculate the size of the movie from the screen size and roi dimensions
% if magCoef is given, use it unless bigger than the maximum allowed by the
% screen size
scrsz = get(0,'ScreenSize');
screenW=scrsz(3);
screenL=scrsz(4);

%if nargin<5 || isempty(magCoef)
magCoef=4;
%end

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


for iFrame=startFrame:endFrame-1
    
    fileNameIm = [char(listOfImages(iFrame,2)) filesep char(listOfImages(iFrame,1))];
    img = double(imread(fileNameIm));

    imagesc(img(minY:maxY,minX:maxX))
    colormap gray

    hold on

    % correct coordinates to correspond to crop from ROI
    xCoord=projData.xCoord(:,iFrame); yCoord=projData.yCoord(:,iFrame);
    outOfRangeIdx=find(xCoord<minX | xCoord>maxX | yCoord<minY | yCoord>maxY);
    xCoord(outOfRangeIdx) = [];
    yCoord(outOfRangeIdx) = [];
    xCoord=xCoord-minX+1;
    yCoord=yCoord-minY+1;

    % get closest colormap index for each feature
    D=createDistanceMatrix(trackVelMatrix(:,iFrame),mapper);
    [sD,idx]=sort(abs(D),2);
    idx(outOfRangeIdx,:)=[];

    scatter(xCoord,yCoord,'Marker','.','cData',cMap(idx(:,1),:));
    
    text(.25,.25,num2str(iFrame),'Color','w','FontWeight','bold','HorizontalAlignment','right','Units','inches')
    
    if aviInstead==1
        F(iFrame) = getframe;
    else
        MakeQTMovie addaxes
        MakeQTMovie('framerate', 5);
    end
end

if aviInstead==1
    movie2avi(F,[movDir filesep movieName '.avi'],'COMPRESSION','Cinepak','FPS',5)
else
    MakeQTMovie finish
end

save([movieName '_roiYX'],'roiYX')
close all
cd(homeDir)
