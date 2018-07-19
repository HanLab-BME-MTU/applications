%This is a script file that displays the colormap of identified adhesion 
% coefficient overlayed on the original cell image.

startFrmNo = imgIndexOfDTimePts(curDTimePt)-firstImgIndex+1;
if length(numAvgFrames) == numDTimePts
   curNumAvgFrames = numAvgFrames(curDTimePt);
else
   curNumAvgFrames = numAvgFrames;
end

%Get the overlaid images of selected image channel.
stackedImg = double(imread(imgFileList{imgChannel}{startFrmNo}));
for k = 1:curNumAvgFrames-1
   stackedImg = stackedImg + double(imread(imgFileList{imgChannel}{k+startFrmNo}));
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

numColors = 128;

%Blue-Green to yellow.
cMap = [linspace(0,1,numColors/2).' ...
   linspace(0.5,1,numColors/2).' ...
   linspace(0.5,0.05,numColors/2).'];
%Yellow to red.
cMap = [cMap; [ones(numColors/2,1) ...
   linspace(1,0,numColors/2).' zeros(numColors/2,1)]];

%Load adhesion coefficients map.
indexStr = sprintf(imgIndexForm,imgIndexOfDTimePts(jj));
adcMapFile = [reslDir filesep 'adcMap' filesep 'adcMap' indexStr '.mat'];
load(adcMapFile);

%Load the saved boundary of the mixed zone.
mixfMapFile = [reslDir filesep 'mixfMap' filesep 'mixfMap' indexStr '.mat'];
load(mixfMapFile);
mixZoneInd = find(~isnan(mixfMap));

if strcmp(showMixZone,'no')
   adcMap(mixZoneInd) = NaN;
end

iDispFieldFileName = ['iDispField' sprintf(imgIndexForm,imgIndex) '.mat'];
iDispFieldFile     = [iDispFieldDir filesep iDispFieldFileName];
s = load(iDispFieldFile);
iDispField = s.iDispField;

gridX = iDispField.gridX;
gridY = iDispField.gridY;

numInd = find(~isnan(adcMap));
maxADC = adfColorDispRange(2)*max(adcMap(numInd));
minADC = adfColorDispRange(1)*min(adcMap(numInd));
adcImg = imDataMapOverlay(stackedImg,adcMap,[0 maxADC],cMap);

figure; hold off;
imshow(adcImg,[]); hold on;

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

