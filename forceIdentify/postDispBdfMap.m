%This is a script file that displays the colormap of identified forces 
% overlayed on the original cell image.

imgIndex = imgIndexOfDTimePts(curDTimePt);

frameNo = imgIndexOfDTimePts(curDTimePt)-firstImgIndex+1;
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

%Load the saved body force map.
indexStr = sprintf(imgIndexForm,imgIndex);
bdfMapFile = [reslDir filesep 'bdfMap' filesep 'bdfMap' indexStr '.mat'];
load(bdfMapFile);

% cMap = colormap('jet');
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

numInd = find(~isnan(bdfMap));
maxBDF = bdfColorDispRange(2)*max(bdfMap(numInd));
minBDF = bdfColorDispRange(1)*max(bdfMap(numInd));
bdfImg = imDataMapOverlay(stackedImg,bdfMap,[minBDF maxBDF],cMap);

figure; hold off;
imshow(bdfImg,[]); hold on;

if strcmp(markMcfBnd,'yes')
   %Load the saved mix zone force map.
   mcfMapFile = [reslDir filesep 'mcfMap' filesep 'mcfMap' indexStr '.mat'];
   load(mcfMapFile);
   mcfZoneInd = find(~isnan(mcfMap));

   %identify the boundary of the mixed zone and draw it.
   mcfBW = zeros(size(mcfMap));
   mcfBW(mcfZoneInd) = 1;
   mcfZoneBnd = bwboundaries(mcfBW);
   for k = 1:length(mcfZoneBnd)
      bndH = plot(mcfZoneBnd{k}(:,2),mcfZoneBnd{k}(:,1));
      set(bndH,'Color',mcfBndColor,'LineWidth',1);
   end
end

if strcmp(markAdfBnd,'yes')
   %Load the saved mix zone force map.
   adfMapFile = [reslDir filesep 'adfMap' filesep 'adfMap' indexStr '.mat'];
   load(adfMapFile);
   adfZoneInd = find(~isnan(adfMap));

   %identify the boundary of the mixed zone and draw it.
   adfBW = zeros(size(adfMap));
   adfBW(adfZoneInd) = 1;
   adfZoneBnd = bwboundaries(adfBW);
   for k = 1:length(adfZoneBnd)
      bndH = plot(adfZoneBnd{k}(:,2),adfZoneBnd{k}(:,1));
      set(bndH,'Color',adfBndColor,'LineWidth',1);
   end
end

if strcmp(markMixZone,'yes')
   %Load the saved mix zone force map.
   mixfMapFile = [reslDir filesep 'mixfMap' filesep 'mixfMap' indexStr '.mat'];
   load(mixfMapFile);
   mixZoneInd = find(~isnan(mixfMap));

   %identify the boundary of the mixed zone and draw it.
   mixBW = zeros(size(mixfMap));
   mixBW(mixZoneInd) = 1;
   mixZoneBnd = bwboundaries(mixBW);
   for k = 1:length(mixZoneBnd)
      bndH = plot(mixZoneBnd{k}(:,2),mixZoneBnd{k}(:,1));
      set(bndH,'Color',mixfBndColor,'LineWidth',1);
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
recDispU1   = forceField.rv(:,1);
recDispU2   = forceField.rv(:,2);
recBFx      = forceField.f(:,1);
recBFy      = forceField.f(:,2);
if strcmp(showFlowVec,'yes')
   quiver(bfDisplayPx,bfDisplayPy, ...
      recDispU1*dispScale,recDispU2*dispScale,0,'y');
end

if strcmp(showBdfVec,'yes')
   quiver(bfDisplayPx,bfDisplayPy, ...
      recBFx*bfScale,recBFy*bfScale,0,'r');
end

iDispFieldFileName = ['iDispField' sprintf(imgIndexForm,imgIndex) '.mat'];
iDispFieldFile     = [iDispFieldDir filesep iDispFieldFileName];
s = load(iDispFieldFile);
iDispField = s.iDispField;

rawDispFieldFileName = ['rawDispField' sprintf(imgIndexForm,imgIndex) '.mat'];
rawDispFieldFile     = [rawDispFieldDir filesep rawDispFieldFileName];
s = load(rawDispFieldFile);
rawDispField = s.rawDispField;


