%This is a script file that displays the colormap of mixed zone of adhesion and myosin forces. 
% overlayed on the original cell image.

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
mixfMapFile = [reslDir filesep 'mixfMap' filesep 'mixfMap' indexStr '.mat'];
load(mixfMapFile);

% cMap = colormap('jet');
cMap = [[zeros(10,1) linspace(0,0.95,10).' linspace(1,0.05,10).']]; %Blue
cMap = [cMap; [ones(30,1) linspace(1,0,30).' zeros(30,1)]]; %Yellow to red.

maxMixF = max(mixfMap(find(~isnan(mixfMap))));
mixImg  = imDataMapOverlay(stackedImg,mixfMap,[0 maxMixF],cMap);
figure; hold off;
imshow(mixImg,[]);

