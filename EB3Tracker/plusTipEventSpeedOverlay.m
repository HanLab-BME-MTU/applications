function plusTipEventSpeedOverlay(projData,velLimit)
% plusTipEventSpeedOverlay creates pause/shrinkage events and growth speed overlays

% INPUT : projData: project data, stored in roi_x/meta folder
%         velLimit: max speed to use for growth speed plot
% OUTPUT: 3 figures:
%         1) pause initiation in yellow and shrinkage initiation in red
%         2) growth speed distribution (so you can tweak velLimit)
%         3) growth speed overlay color-coded by speed

if nargin<1
    error('plusTipEventSpeedOverlay: not enough input arguments')
end
if nargin<2 || isempty(velLimit)
    velLimit=30;
end

% format image/analysis directory paths
imDir=formatPath(projData.imDir);
anDir=formatPath(projData.anDir);

% get first image from imDir
[listOfImages] = searchFiles('.tif',[],imDir,0);
img = double(imread([listOfImages{1,2} filesep listOfImages{1,1}]));

% extract track info (aggregated)
allData=abs(projData.nTrack_sF_eF_vMicPerMin_trackType_lifetime_totalDispPix);

% pause info
pauseIdx=find(allData(:,5)==2);  % index of pause events
growthRow=allData(pauseIdx,1);   % track ID
growthCol=allData(pauseIdx+1,2); % frame where growth phase after pause begins
indPause=sub2ind(size(projData.xCoord),growthRow, growthCol);

% shrink info
shrinkIdx=find(allData(:,5)==3);  % index of shrinkage events
growthRow=allData(shrinkIdx,1);   % track ID
growthCol=allData(shrinkIdx+1,2); % frame where growth phase after shrink begins
indShrink=sub2ind(size(projData.xCoord),growthRow, growthCol);


% get size of the images
minY=1; maxY=size(img,1);
minX=1; maxX=size(img,2);
scrsz = get(0,'ScreenSize');
screenW=scrsz(3);
screenL=scrsz(4);
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
figPos=[round(screenW*(1-movieW/screenW)/2) round(screenL*(1-movieL/screenL)/2) movieW movieL];




% plot pause in yellow and shrinkage in red
h=figure('Position',figPos);
imagesc(img); colormap gray;
hold on
scatter(projData.xCoord(indPause),projData.yCoord(indPause),'y','filled')
scatter(projData.xCoord(indShrink),projData.yCoord(indShrink),'r','filled')
xlabel('pause (yellow) and shrinkage (red) initiation sites')


% extract growth phases only
trackType=allData(:,5);
[xMat,yMat]=plusTipGetSubtrackCoords(projData,[]);
xMat=xMat(trackType==1,:);
yMat=yMat(trackType==1,:);
vel=allData(trackType==1,4);

% check upper limit of growth speed
figure; hist(vel,25)
title('growth speed distribution (um/min)');


cMapLength=128; cMap=jet(cMapLength);
m=max(abs([min(vel); max(vel)]));
m=round(min(m,velLimit));
vel(vel<-m)=-m;
vel(vel>m)=m;
mapper=linspace(0,m,cMapLength)';


x1=xMat(:,1:end-1)'; 
y1=yMat(:,1:end-1)'; 


% get closest colormap index for each feature
D=createDistanceMatrix(vel,mapper);
[sD,idx]=sort(abs(D),2);

h=figure('Position',figPos);
imagesc(img); colormap gray;
hold on
xlabel(['speed map of growth trajectories, max velocity ' num2str(velLimit) ' um/min'])
for i=1:size(x1,2)
    plot(x1(:,i),y1(:,i),'color',cMap(idx((i),1),:));
end




