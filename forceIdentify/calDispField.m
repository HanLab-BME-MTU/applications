%Preprocessing of the raw experimental data. It creates by interpolation of
% the raw data a displacement field on grid points inside a predefined polygon.

%Load parameters that are needed for the calculation of the displacement
% field.
run([resultPath 'setPar']);

%Get the image of the cell.
cellImg = imread(imgFile{1});

%Show the image 
figure(gcf); hold off;
imshow(cellImg,[]); axis on; hold on;

%Load the displacement field given in the multidimentional array, 'M' that is
% produced by the speckle tracking software.
load(dataFile);

%Load the polygon, 'fieldPG' that define the geometry of the lamollipodium.
load([modelPath 'fieldGeom']);

%Plot the polygon of the lamollipodium.
plot(fieldPGx,fieldPGy,'r-');

rawDispV  = cell(numTimeSteps,1);
rawDataU1 = cell(numTimeSteps,1);
rawDataU2 = cell(numTimeSteps,1);

dataPx   = cell(numTimeSteps,1);
dataPy   = cell(numTimeSteps,1);
speckleP = cell(numTimeSteps,1);

%Generate a grid with gridDx X gridDy raster sidelength. And, 
% the grid sampling data points are those inside the predefined polygon.
[numPixelsY,numPixelsX] = size(cellImg);
numGridX = ceil(numPixelsX/gridDx);
numGridY = ceil(numPixelsY/gridDy);
gridx    = linspace(1,numPixelsX,numGridX);
gridy    = linspace(1,numPixelsY,numGridY);
[gridX,gridY] = meshgrid(gridx,gridy);
gridX = reshape(gridX,length(gridX(:)),1);
gridY = reshape(gridY,length(gridY(:)),1);
in   = inpolygon(gridX,gridY,fieldPGx,fieldPGy);
ind  = find(in==1);

gridPx = gridX(ind);
gridPy = gridY(ind);

frame1 = startFrame;
frame2 = frame1+framesPerTraj-1;
for jj = 1:numTimeSteps
   speckleP{jj} = M(find(M(:,1,frame1)~=0 & M(:,3,frame1)~=0),1:2,frame1);
   if strcmp(dataSite,'grid')
      dataPy{jj} = gridPy;
      dataPx{jj} = gridPx;
   elseif strcmp(dataSite,'speckle')
      in    = inpolygon(speckleP{jj}(:,2),speckleP{jj}(:,1),fieldPGx,fieldPGy);
      ind   = find(in==1);

      dataPy{jj} = speckleP{jj}(ind,1);
      dataPx{jj} = speckleP{jj}(ind,2);
   end

   if strcmp(dispType,'MFAverage')
      rawDispV{jj}  = M(find(M(:,1,frame1)~=0 & M(:,3,frame1)~=0),:,frame1);
      %rawDataPx{jj} = rawDispV{jj}(:,2);
      %rawDataPy{jj} = rawDispV{jj}(:,1);
      for k = frame1+1:frame2-1
         rawDispV{jj} = [rawDispV{jj}; ...
            M(find(M(:,1,k)~=0 & M(:,3,k)~=0),:,k)];
      end
      dataDispV = vectorFieldSparseInterp(rawDispV{jj}, ...
         [dataPy{jj} dataPx{jj}],2*MFCorLen,MFCorLen,[]); 
   elseif strcmp(dispType,'MFTrack')
      %MFT          = mFrameTrajBuild(M(:,:,frame1:frame2-1));
      MFT = mFrameTrajBuild(M(:,:,frame1:frame2-1),MFCorLen, ...
         [dataPy{jj} dataPx{jj}]);
      dataDispV    = MFT(:,[1 2 end-1 end]);
      rawDispV{jj} = dataDispV;
   elseif strcmp(dispType,'SFrame')
      rawDispV{jj} = M(find(M(:,1,frame1)~=0 & M(:,3,frame1)~=0),:,frame1);
      if strcmp(dataSite,'grid')
         dataDispV = vectorFieldSparseInterp(rawDispV{jj}, ...
            [dataPy{jj} dataPx{jj}],2*MFCorLen,MFCorLen,[]); 
      else
         dataDispV = rawDispV{jj};
      end
   end
   frame1 = frame1+framesPerTStep;
   frame2 = frame1+framesPerTraj-1;

   %Extract data points and the raw displamenents on the data points.
   dataPy{jj} = dataDispV(:,1);
   dataPx{jj} = dataDispV(:,2);
   rawDataU1{jj} = dataDispV(:,4)-dataDispV(:,2);
   rawDataU2{jj} = dataDispV(:,3)-dataDispV(:,1);
end

%Plot the data points.
%plot(dataPx{1},dataPy{1},'y.'); 

%Plot the data points and displacements.
quiver(dataPx{1},dataPy{1},rawDataU1{1}*dispScale, ...
   rawDataU2{1}*dispScale,0,'y'); 

sDataU1 = cell(numTimeSteps,1);
sDataU2 = cell(numTimeSteps,1);

if strcmp(calInterp,'yes')
   for jj = 1:numTimeSteps
      %Smoothed displacements at the data points.
      sDispV = vectorFieldSparseInterp(rawDispV{jj}, ...
         [dataPy{jj} dataPx{jj}], 2*corLen,corLen,[fieldPGx fieldPGy]); 
      sDataU1{jj} = sDispV(:,4)-sDispV(:,2);
      sDataU2{jj} = sDispV(:,3)-sDispV(:,1);
   end

   if strcmp(showInterp,'yes')
      %the filtered displacements.
      quiver(dataPx{1},dataPy{1},sDataU1{1}*dispScale, ...
         sDataU2{1}*dispScale,0,'r'); 
   end

   save([resultPath 'dispField'],'rawDispV','dataPx', ...
      'dataPy','rawDataU1','rawDataU2','sDataU1','sDataU2', ...
      'MFCorLen','corLen','edgCorLen','gridPx','gridPy', ...
      'gridx','gridy','gridX','gridY','speckleP');
else
   save([resultPath 'dispField'],'rawDispV','dataPx', ...
      'dataPy','rawDataU1','rawDataU2','MFCorLen','corLen','edgCorLen', ...
      'gridPx','gridPy','gridx','gridy','gridX','gridY','speckleP');
end


hold off;
