%This is a script file that displays the colormap of the flow speed
% overlaid to the original cell image.

startImgIndex = imgIndexOfDTimePts(curDTimePt);
startFrmNo    = startImgIndex-firstImgIndex+1;
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

%Load the saved myosin contraction force map.
indexStr = sprintf(indexForm,imgIndexOfDTimePts(jj));
spdMapFile = [reslDir filesep 'spdMap' filesep 'spdMap' indexStr '.mat'];
load(spdMapFile);

%cMap = colormap('jet');
%cMap = [[zeros(10,1) linspace(0,0.95,10).' linspace(1,0.05,10).']]; %Blue
%cMap = [cMap; [ones(30,1) linspace(1,0,30).' zeros(30,1)]]; %Yellow to red.
%Blue-Green to yellow.
numSpdColors = 128;
cMap = [linspace(0,1,numSpdColors/2).' ...
   linspace(0.5,1,numSpdColors/2).' ...
   linspace(0.5,0.05,numSpdColors/2).'];
%Yellow to red.
cMap = [cMap; [ones(numSpdColors/2,1) ...
   linspace(1,0,numSpdColors/2).' zeros(numSpdColors/2,1)]];


spdMap = spdMap*spdUnitConvFactor;
maxSpd = spdColorDispRange(2)*max(spdMap(find(~isnan(spdMap))));
minSpd = spdColorDispRange(1)*max(spdMap(find(~isnan(spdMap))));
spdMap(find(spdMap>maxSpd)) = maxSpd;
spdMap(find(spdMap<minSpd)) = minSpd;

figure; hold off;
if strcmp(showSpdCBar,'yes')
   mapH = imshow(spdMap,[]); colormap(cMap); colorbar; hold on;
   delete(mapH);
end
spdImg = imDataMapOverlay(stackedImg,spdMap,[minSpd maxSpd],cMap);
imshow(spdImg,[]);

if strcmp(showFlowVec,'yes')
   jj = curDTimePt;
   quiver(dataPx{jj},dataPy{jj}, ...
      dataUC1{jj}*dispScale,dataUC2{jj}*dispScale,0,'y');
end

