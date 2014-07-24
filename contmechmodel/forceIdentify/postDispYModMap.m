%This is a script file that displays the colormap of the Young's modulus used in
% the force reconstruction. It is overlaid to the original cell image.

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

%Get the dimension of the cell image.
pixelX = [1:size(stackedImg,2)];
pixelY = [1:size(stackedImg,1)];

%Load the saved myosin contraction force map.
indexStr = sprintf(indexForm,imgIndexOfDTimePts(jj));
ymodMapFile = [reslDir filesep 'ymodMap' filesep 'ymodMap' indexStr '.mat'];
load(ymodMapFile);

cMap = colormap('jet');
%cMap = [[zeros(10,1) linspace(0,0.95,10).' linspace(1,0.05,10).']]; %Blue
%cMap = [cMap; [ones(30,1) linspace(1,0,30).' zeros(30,1)]]; %Yellow to red.

figure; hold off;

if strcmp(showYModCBar,'yes')
   mapH = imshow(ymodMap,[]);
   colormap(cMap); 
   colorbar; hold on;
   delete(mapH);
end

maxYMod = max(ymodMap(find(~isnan(ymodMap))));
ymodImg = imDataMapOverlay(stackedImg,ymodMap,[0 maxYMod],cMap,1);
imshow(ymodImg,[]);

