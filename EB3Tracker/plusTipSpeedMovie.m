function plusTipSpeedMovie(projData,timeRange,velLimit,roiYX,aviInstead)
% plusTipSpeedMovie makes a movie of features color-coded by velocity
%
%
%OUTPUT Quicktime movies and the regions of interest used to generate them.
%       
%INPUT  projData          : output of metaEB3analysis, stored in /meta
%                           if not given, user will be asked to select the
%                           file
%       timeRange         : row vector of the form [startFrame endFrame]
%                           indicating time range to plot. if not given or
%                           given as [], tracks from the whole movie will
%                           be displayed. if indivTrack is given, this
%                           parameter defaults to the frames over which the
%                           track appears, regardless of user input.
%       velLimit          : max speed to use for colormap (speeds higher
%                           than velLimit or lower than -velLimit will be
%                           plotted as red or blue, respectively. if not
%                           given or given as [], max value will be used.
%       roiYX             : coordinates of a region-of-interest (closed
%                           polygon) in which tracks should be plotted, or
%                           a logical mask the size of the image. if
%                           indivTrack is given, this parameter
%                           defaults to the ROI circumscribing the full
%                           track, regardless of user input. if given as
%                           [], user will be asked to select a ROI.
%       aviInstead        : 1  to make AVI movie (works in Windows only),
%                          {0} to make MOV movie
%
%OUTPUT One or more movies and the regions of interest used to
%       generate them. Circle are growths, triangles are fgaps, and squares
%       are bgaps

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
if ~isfield(projData,'movDir')
    projData.movDir=uigetdir(projData.anDir,'Please select OUTPUT directory');
end
if isequal(projData.movDir,0)
    disp('No output directory selected.')
    return
end
cd(projData.movDir)

% check whether a time range for plotting was input
if nargin<3 || isempty(timeRange)
    timeRange = [1 projData.numFrames];
else
    if timeRange(1) < 1
        timeRange(1)=1;
    end
    if timeRange(2) > projData.numFrames
        timeRange(2)=projData.numFrames-1;
    end
end
startFrame = timeRange(1);
endFrame   = timeRange(2);
    

% default: use maximum velocity as limit
if nargin<3 || isempty(velLimit) || isinf(velLimit)
    velLimit = [];
end

% get list of files and image size
[listOfImages] = searchFiles('.tif',[],projData.imDir,0);
fileNameIm = [char(listOfImages(1,2)) filesep char(listOfImages(1,1))];
img = double(imread(fileNameIm));
[imL,imW] = size(img);


% check input for ROI coordinates
if nargin<4 || isempty(roiYX)
    imscaled=(img-min(img(:)))./(max(img(:))-min(img(:)));
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

% use abs of instantaneous velocities, since bgaps have negative values
trackVelMatrix=abs(projData.frame2frameVel_micPerMin);

% get merged matrix so all we have are growths, fgaps, and bgaps (1,2,3)
dataMatMerge=plusTipMergeSubtracks(projData); % first output is merged

% make corresponding matrix to velocity matrix where values represent type
trackTypeMtrix=zeros(size(projData.xCoord));
for iSub=size(dataMatMerge,1):-1:1
    sF=dataMatMerge(iSub,2);
    eF=dataMatMerge(iSub,3);
    fullIdx=dataMatMerge(iSub,1);
    trackType=dataMatMerge(iSub,5);
    trackTypeMtrix(fullIdx,sF:eF)=trackType;
end

% get velocity limit for color map
cMapLength=128; 
cMap=mat2cell(jet(cMapLength),ones(cMapLength,1),3);

% find max velocity value m
if ~isempty(velLimit)
    m=round(velLimit);
else
    m=round(nanmax(nanmax(trackVelMatrix(:,startFrame:endFrame))));
end
trackVelMatrix(trackVelMatrix>m)=m;
mapper=linspace(0,m,cMapLength)';


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
while exist([projData.movDir filesep temp '.mov'],'file')>0
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

% magnification coefficient
magCoef=inf;

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

if aviInstead==0
    eval(['MakeQTMovie start ', movieName '.mov']);
end

frmCounter=1;
mrkTpe={'o';'^';'s'};
prop_name(1) = {'Marker'};
prop_name(2) = {'MarkerFaceColor'};
prop_name(3) = {'MarkerEdgeColor'};
for iFrame=startFrame:endFrame-1
    
    fileNameIm = [char(listOfImages(iFrame,2)) filesep char(listOfImages(iFrame,1))];
    img = double(imread(fileNameIm));
    imagesc(img(minY:maxY,minX:maxX))
    colormap gray

    hold on

    % correct coordinates to correspond to crop from ROI
    xCoord=projData.xCoord(:,iFrame); 
    yCoord=projData.yCoord(:,iFrame);
    trckTp=trackTypeMtrix(:,iFrame);
    outOfRangeIdx=find(xCoord<minX | xCoord>maxX | yCoord<minY | yCoord>maxY | isnan(xCoord));
    xCoord(outOfRangeIdx) = [];
    yCoord(outOfRangeIdx) = [];
    trckTp(outOfRangeIdx) = [];
    xCoord=xCoord-minX+1;
    yCoord=yCoord-minY+1;

    
    % get closest colormap index for each feature
    vel=trackVelMatrix(:,iFrame);
    vel(outOfRangeIdx) = [];
    D=createDistanceMatrix(vel,mapper);
    [sD,idx]=sort(abs(D),2);
    idx=idx(:,1);
    
    % record marker type and face/edge colors for each feature
    nC=length(xCoord);
    prop_values(1:nC,1) = mrkTpe(trckTp);
    prop_values(1:nC,2) = cMap(idx,:);
    prop_values(1:nC,3) = cMap(idx,:);
    
    % use plot instead of scatter so more flexibility with properties. to
    % do this, make 2 x nPoints matrix where the second row is all NaNs and
    % then use the plot function
    xCoord=[xCoord nan(size(xCoord))]';
    yCoord=[yCoord nan(size(yCoord))]';
    h=plot(xCoord,yCoord);
    set(h,prop_name,prop_values)

    % these have to be defined each frame
    clear prop_values

    text(.25,.25,num2str(iFrame),'Color','w','FontWeight','bold','HorizontalAlignment','right','Units','inches')
    
    if aviInstead==1
        F(frmCounter) = getframe;
    else
        MakeQTMovie addaxes
        MakeQTMovie('framerate', 5);
    end
    frmCounter=frmCounter+1;
end

if aviInstead==1
    movie2aviNADA_CAW(F,[projData.movDir filesep movieName '.avi'],'COMPRESSION','Cinepak','FPS',5)
else
    MakeQTMovie finish
end

save([movieName '_roiYX'],'roiYX')

cd(homeDir)
close(figure(gcf))
