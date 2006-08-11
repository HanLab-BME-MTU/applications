%This is a script file that displays the colormap of identified myosin contraction
% overlayed on the original cell image.

if ~isdir([reslDir filesep 'mcfFig'])
   success = mkdir(reslDir,'mcfFig');
   if ~success
      error('trouble making directory.')
   end
end
mcfFigDir = [reslDir filesep 'mcfFig'];

mcfTifDir = [reslDir filesep 'mcfTif'];
if ~isdir([reslDir filesep 'mcfTif'])
   success = mkdir(reslDir,'mcfTif');
   if ~success
      error('trouble making directory.')
   end
end
mcfTifDir = [reslDir filesep 'mcfTif'];

answer = input('Select time steps (0 for all):');
if isempty(answer) | answer == 0
   selTimeSteps = 1:numDTimePts;
else
   selTimeSteps = answer;
end

%Find the maximum force.
maxF = 0;
for ii = 1:length(selTimeSteps)
   jj = selTimeSteps(ii);
   
   if jj == -1
      imgIndex = [];
      dispImgIndex = firstImgIndex;
   else
      imgIndex = imgIndexOfDTimePts(jj);
      dispImgIndex = max(firstImgIndex,imgIndex);
   end
   
   %Load the saved body force map.
   indexStr = sprintf(imgIndexForm,imgIndex);
   bdfMapFile = [reslDir filesep 'bdfMap' filesep 'bdfMap' indexStr '.mat'];
   load(bdfMapFile);
   
   numInd = find(~isnan(bdfMap));
   maxF = max(maxF,max(bdfMap(numInd)));
end

%Creat color map
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

figH = figure;
backStr = '';
for ii = 1:length(selTimeSteps)
   jj = selTimeSteps(ii);

   localStartTime = cputime;

   for kk = 1:length(backStr)
      fprintf(1,'\b');
   end
   backStr = sprintf('   Time step %d ... ',jj);
   fprintf(1,backStr);

   if jj == -1
      imgIndex = [];
      dispImgIndex = firstImgIndex;
   else
      imgIndex = imgIndexOfDTimePts(jj);
      dispImgIndex = max(firstImgIndex,imgIndex);
   end

   frameNo = dispImgIndex-firstImgIndex+1;
   if length(numAvgFrames) == numDTimePts
      if jj == -1
         curNumAvgFrames = numAvgFrames(1);
      else
         curNumAvgFrames = numAvgFrames(jj);
      end
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

   %Load the saved myosin contraction force map.
   indexStr = sprintf(imgIndexForm,imgIndex);
   mcfMapFile = [reslDir filesep 'mcfMap' filesep 'mcfMap' indexStr '.mat'];
   load(mcfMapFile);

   %Load the saved mix zone force map.
   mixfMapFile = [reslDir filesep 'mixfMap' filesep 'mixfMap' indexStr '.mat'];
   load(mixfMapFile);
   mixZoneInd = find(~isnan(mixfMap));

   if strcmp(showMixZone,'no')
      mcfMap(mixZoneInd) = NaN;
   end

   iDispFieldFileName = ['iDispField' sprintf(imgIndexForm,imgIndex) '.mat'];
   iDispFieldFile     = [iDispFieldDir filesep iDispFieldFileName];
   s = load(iDispFieldFile);
   iDispField = s.iDispField;

   gridX = iDispField.gridX;
   gridY = iDispField.gridY;

   maxMCF = mcfColorDispRange(2)*maxF;
   minMCF = mcfColorDispRange(1)*maxF;
   mcfImg = imDataMapOverlay(stackedImg,mcfMap,[minMCF maxMCF],cMap);

   figure(figH); hold off;
   imshow(mcfImg,[]); hold on;

   if strcmp(markMixZone,'yes')
      %identify the boundary of the mixed zone and draw it.
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

   recDispU1 = iDispField.rv(:,1);
   recDispU2 = iDispField.rv(:,2);

   if strcmp(showFlowVec,'yes')
      quiver(bfDisplayPx,bfDisplayPy, ...
         recDispU1*dispScale,recDispU2*dispScale,0,'y');
   end

   if strcmp(showMcfVec,'yes')
      bfDisplayPxPix = round(bfDisplayPx);
      bfDisplayPyPix = round(bfDisplayPy);
      [m,n] = size(mcfMap);
      bfDisplayPxPix(find(bfDisplayPxPix>n)) = n;
      bfDisplayPxPix(find(bfDisplayPxPix<1)) = 1;
      bfDisplayPyPix(find(bfDisplayPyPix>m)) = m;
      bfDisplayPyPix(find(bfDisplayPyPix<1)) = 1;
      dispInd = sub2ind(size(mcfMap),bfDisplayPyPix,bfDisplayPxPix);
      mcfVecInd = find(~isnan(mcfMap(dispInd)) & isnan(mixfMap(dispInd)));
      quiver(bfDisplayPx(mcfVecInd),bfDisplayPy(mcfVecInd), ...
         recBFx(mcfVecInd)*bfScale,recBFy(mcfVecInd)*bfScale,0,'r');
   end
   titleStr = sprintf(['Contraction Force\n' 'Time Step: %d Image Index: %d'], ...
      jj,imgIndex);
   title(titleStr);
   
   %Save the figure
   mcfFigFile = [mcfFigDir filesep 'mcfFig' indexStr '.fig'];
   saveas(figH,mcfFigFile,'fig');
   mcfTifFile = [mcfTifDir filesep 'mcfTif' indexStr '.tif'];
   print(figH,'-dtiffnocompression',mcfTifFile);   
end
fprintf(1,'\n');

