%This is a script file that displays the colormap of identified adhesion sites
% overlayed on the original cell image.

imgIndex = imgIndexOfDTimePts(curDTimePt);

frameNo = imgIndex-firstImgIndex+1;
if length(numAvgFrames) == numDTimePts
   curNumAvgFrames = numAvgFrames(curDTimePt);
else
   curNumAvgFrames = numAvgFrames;
end

%Get the overlaid images of selected image channel.
stackedImg = double(imread(imgFileList{imgChannel}{frameNo}));
for k = 1:curNumAvgFrames-1
   stackedImg = stackedImg + double(imread(imgFileList{imgChannel}{k+frameNo}));
end
stackedImg = stackedImg/curNumAvgFrames;

maxImgI = max(stackedImg(:));
minImgI = min(stackedImg(:));
if minImgI < maxImgI
   stackedImg = (stackedImg-minImgI)/(maxImgI-minImgI);
   stackedImg(find(stackedImg>imgIRange(2))) = imgIRange(2);
   stackedImg(find(stackedImg<imgIRange(1))) = imgIRange(1);
end

%Get the dimension of the cell image.
pixelX = [1:size(stackedImg,2)];
pixelY = [1:size(stackedImg,1)];

%Load adhesion force.
indexStr = sprintf(imgIndexForm,imgIndex);
adfMapFile = [reslDir filesep 'adfMap' filesep 'adfMap' indexStr '.mat'];
load(adfMapFile);

%Load the boundary of the mixed zone.
mixfMapFile = [reslDir filesep 'mixfMap' filesep 'mixfMap' indexStr '.mat'];
load(mixfMapFile);
mixZoneInd = find(~isnan(mixfMap));

if strcmp(showMixZone,'no')
   adfMap(mixZoneInd) = NaN;
end

%cMap = [[zeros(10,1) linspace(0,0.95,10).' linspace(1,0.05,10).']]; %Blue
%cMap = [cMap; [ones(30,1) linspace(1,0,30).' zeros(30,1)]]; %Yellow to red.
numColors = 128;

%Blue-Green to yellow.
cMap = [linspace(0,1,numColors/2).' ...
   linspace(0.5,1,numColors/2).' ...
   linspace(0.5,0.05,numColors/2).'];
%Yellow to red.
cMap = [cMap; [ones(numColors/2,1) ...
   linspace(1,0,numColors/2).' zeros(numColors/2,1)]];

numInd = find(~isnan(adfMap));
maxADF = adfColorDispRange(2)*max(adfMap(numInd));
minADF = adfColorDispRange(1)*max(adfMap(numInd));
adfImg = imDataMapOverlay(stackedImg,adfMap,[minADF maxADF],cMap);

figure; hold off;
imshow(adfImg,[]); hold on;

if strcmp(markMixZone,'yes')
   mixBW = zeros(size(mixfMap));
   mixBW(mixZoneInd) = 1;
   mixZoneBnd = bwboundaries(mixBW);
   for k = 1:length(mixZoneBnd)
      plot(mixZoneBnd{k}(:,2),mixZoneBnd{k}(:,1),'w','LineWidth',1);
   end

   %Also draw some grey dots in the mixed region.
   gridInd = sub2ind(size(mixBW),round(gridY),round(gridX));
   mixGridInd = find(mixBW(gridInd)==1);
   dotH = plot(gridX(mixGridInd),gridY(mixGridInd),'w.');
   set(dotH,'Color',[0.7 0.7 0.7],'LineWidth',3);
end

forceFieldFile = [forceFieldDir filesep 'forceField' ...
   sprintf(imgIndexForm,imgIndex) '.mat'];
s = load(forceFieldFile);
forceField = s.forceField;

bfDisplayPx = forceField.p(:,1);
bfDisplayPy = forceField.p(:,2);
recBFx      = forceField.f(:,1);
recBFy      = forceField.f(:,2);

iDispFieldFileName = ['iDispField' sprintf(imgIndexForm,imgIndex) '.mat'];
iDispFieldFile     = [iDispFieldDir filesep iDispFieldFileName];
s = load(iDispFieldFile);
iDispField = s.iDispField;

recDispU1   = iDispField.rv(:,1);
recDispU2   = iDispField.rv(:,2);
if strcmp(showFlowVec,'yes')
   quiver(bfDisplayPx,bfDisplayPy, ...
      recDispU1*dispScale,recDispU2*dispScale,0,'y');
end

if strcmp(showAdfVec,'yes')
   bfDisplayPxPix = round(bfDisplayPx);
   bfDisplayPyPix = round(bfDisplayPy);
   [m,n] = size(adfMap);
   bfDisplayPxPix(find(bfDisplayPxPix>n)) = n;
   bfDisplayPxPix(find(bfDisplayPxPix<1)) = 1;
   bfDisplayPyPix(find(bfDisplayPyPix>m)) = m;
   bfDisplayPyPix(find(bfDisplayPyPix<1)) = 1;
   dispInd = sub2ind(size(adfMap),bfDisplayPyPix,bfDisplayPxPix);
   adfVecInd = find(~isnan(adfMap(dispInd)));
   quiver(bfDisplayPx(adfVecInd),bfDisplayPy(adfVecInd), ...
      recBFx(adfVecInd)*bfScale,recBFy(adfVecInd)*bfScale,0,'r');
end

