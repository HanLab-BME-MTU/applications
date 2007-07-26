%This is a script file that displays the colormap of the flow speed
% overlaid to the original cell image.

if strcmp(spdUnit,'pixelPerFrame')
   spdUnitConvFactor = 1;
   spdUnitStr        = 'pixel/frame';
elseif strcmp(spdUnit,'pixelPerMin')
   spdUnitConvFactor = 60/frameInterval;
   spdUnitStr        = 'pixel/min';
elseif strcmp(spdUnit,'umPerMin')
   spdUnitConvFactor = 60/frameInterval*actPixelSize/1000;
   spdUnitStr        = 'um/min';
elseif strcmp(spdUnit,'nmPerSec')
   spdUnitConvFactor = 1/frameInterval*actPixelSize;
   spdUnitStr        = 'nm/sec';
end

%Find the maximum speed of the whole stack.
globalMaxSpeed = 0;

if ~isdir([reslDir filesep 'spdFig'])
   success = mkdir(reslDir,'spdFig');
   if ~success
      error('trouble making directory.')
   end
end
spdFigDir = [reslDir filesep 'spdFig'];

if ~isdir([reslDir filesep 'spdTif'])
   success = mkdir(reslDir,'spdTif');
   if ~success
      error('trouble making directory.')
   end
end
spdTifDir = [reslDir filesep 'spdTif'];

answer = input('Select time steps (0 for all):');
if isempty(answer) | answer == 0
   selTimeSteps = 1:numDTimePts;
else
   selTimeSteps = answer;
end

%Find the maximum speed for each of the selected time steps and the overall
% max speed..
maxSpeed = zeros(1,length(selTimeSteps));
for ii = 1:length(selTimeSteps)
   jj = selTimeSteps(ii);
   
   imgIndex = imgIndexOfDTimePts(jj);
   dispImgIndex = max(firstImgIndex,imgIndex);

   %Load the saved body force map.
   indexStr = sprintf(imgIndexForm,imgIndex);
   spdMapFile = [reslDir filesep 'spdMap' filesep 'spdMap' indexStr '.mat'];
   s = load(spdMapFile);
   spdMap = s.spdMap;
   
   numInd = find(~isnan(spdMap));
   maxSpeed(ii) = max(spdMap(numInd));
end
globalMaxSpeed = max(maxSpeed);

if isinf(maxSPDToShow)
   maxDispSPD = spdColorDispRange(2)*globalMaxSpeed;
else
   maxDispSPD = maxSPDToShow;
end
if isinf(minSPDToShow)
   minDispSPD = spdColorDispRange(1)*globalMaxSpeed;
else
   minDispSPD = minSPDToShow;
end

maxDispSPD     = maxDispSPD*spdUnitConvFactor;
minDispSPD     = minDispSPD*spdUnitConvFactor;
maxSpeed       = maxSpeed*spdUnitConvFactor;
globalMaxSpeed = globalMaxSpeed*spdUnitConvFactor;

%Load detected cell edge files;
edge_sp_array_x = [];
edge_sp_array_y = [];
pixel_edge      = [];
if strcmp(markCellEdge,'yes') && isdir(edgeDir)
   edge_splineFile = [edgeDir filesep 'edge_spline.mat'];
   pixel_edgeFile  = [edgeDir filesep 'pixel_edge.mat'];
   if exist(edge_splineFile,'file') && exist(pixel_edgeFile,'file')
      s = load(edge_splineFile);
      edge_sp_array_x = s.edge_sp_array_x;
      edge_sp_array_y = s.edge_sp_array_y;

      s = load(pixel_edgeFile);
      pixel_edge = s.pixel_edge;
   end
end

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

if strcmp(showSpdCBar,'yes')
   figH1 = figure;
end
figH2 = figure;
backStr = '';
for ii = 1:length(selTimeSteps)
   jj = selTimeSteps(ii);

   localStartTime = cputime;

   for kk = 1:length(backStr)
      fprintf(1,'\b');
   end
   backStr = sprintf('   Time step %d ... ',jj);
   fprintf(1,backStr);

   if jj == 0 
      imgIndex = [];
      dispImgIndex = firstImgIndex;
   else
      imgIndex = imgIndexOfDTimePts(jj);
      dispImgIndex = max(firstImgIndex,imgIndex);
   end

   frameNo = dispImgIndex-firstImgIndex+1;
   if length(numAvgFrames) == numDTimePts
      if jj == 0
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

   %Load the saved body force map.
   indexStr = sprintf(imgIndexForm,imgIndex);
   spdMapFile = [reslDir filesep 'spdMap' filesep 'spdMap' indexStr '.mat'];
   load(spdMapFile);
   spdMap   = spdMap*spdUnitConvFactor;

   %Get the dimension of the cell image.
   pixelX = [1:size(stackedImg,2)];
   pixelY = [1:size(stackedImg,1)];

   iDispFieldFileName = ['iDispField' sprintf(imgIndexForm,imgIndex) '.mat'];
   iDispFieldFile     = [iDispFieldDir filesep iDispFieldFileName];
   s = load(iDispFieldFile);
   iDispField = s.iDispField;

   gridX = iDispField.gridX;
   gridY = iDispField.gridY;
  
   %Create color bar figure.
   if strcmp(showSpdCBar,'yes')
      cSpdMap = spdMap;
      cSpdMap(1,1) = globalMaxSpeed;
      cSpdMap(1,2) = 0;
      cSpdMap(find(cSpdMap>=maxDispSPD)) = maxDispSPD;
      cSpdMap(find(cSpdMap<=minDispSPD)) = minDispSPD;
      figure(figH1);
      imshow(cSpdMap,[]); colormap(cMap); colorbar; hold on;
      spdCBarFile = [spdFigDir filesep 'spdCBar' sprintf(imgIndexForm,imgIndex) '.fig'];
      saveas(figH1,spdCBarFile,'fig');
   end

   spdImg = imDataMapOverlay(stackedImg,spdMap,[minDispSPD maxDispSPD],cMap);

   figure(figH2); hold off;
   imshow(spdImg,[]); hold on;

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

   switch bfDisplaySite
      case 'dataSite'
         bfDisplayPx = forceField.p(:,1);
         bfDisplayPy = forceField.p(:,2);
         recBFx      = forceField.f(:,1);
         recBFy      = forceField.f(:,2);

         recDispU1 = iDispField.rv(:,1);
         recDispU2 = iDispField.rv(:,2);
      case 'everyOther'
         bfDisplayPx = forceField.p(1:2:end,1);
         bfDisplayPy = forceField.p(1:2:end,2);
         recBFx      = forceField.f(1:2:end,1);
         recBFy      = forceField.f(1:2:end,2);

         recDispU1 = iDispField.rv(1:2:end,1);
         recDispU2 = iDispField.rv(1:2:end,2);
   end

   if strcmp(showFlowVec,'yes')
      quiver('v6',bfDisplayPx,bfDisplayPy, ...
         recDispU1*dispScale,recDispU2*dispScale,0,'y');
   end

   if strcmp(showBdfVec,'yes')
      quiver('v6',bfDisplayPx,bfDisplayPy, ...
         recBFx*bfScale,recBFy*bfScale,0,'r');
   end

   titleStr = sprintf(['Speed Map\n' 'Time Step: %d Image Index: %d ' ...
      'Residue: %5.2f Max Speed: %f (Unit: ' spdUnit ')'], ...
      jj,imgIndex,iDispField.vRelRes,maxSpeed(ii));
   title(titleStr);
   
   %Save the figure
   spdFigFile = [spdFigDir filesep 'spdFig' indexStr '.fig'];
   saveas(figH2,spdFigFile,'fig');
   
   spdTifFile = [spdTifDir filesep 'spdTif' indexStr '.tif'];
   print(figH2,'-dtiffnocompression',spdTifFile);
end
fprintf(1,'\n');


