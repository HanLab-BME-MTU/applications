%Preprocessing of the raw experimental data. It creates by interpolation of
% the raw data a displacement field on grid points inside a predefined polygon.

%Load parameters that are needed for the calculation of the displacement
% field.
run([resultPath 'setPar']);

%Get the image of the cell.
rawI = cell(numTimeSteps,1);
for jj = 1:numTimeSteps
   rawI{jj} = imread(imgFile{jj});
end

%Show the image 
figure(gcf); hold off;
imshow(rawI{1},[]); axis on; hold on;

%Load the displacement field given in the multidimentional array, 'M' that is
% produced by the speckle tracking software.
load(dataFile);

rawDispV  = cell(numTimeSteps,1);
rawDataPx = cell(numTimeSteps,1);
rawDataPy = cell(numTimeSteps,1);
rawDataU1 = cell(numTimeSteps,1);
rawDataU2 = cell(numTimeSteps,1);

frame1 = startFrame;
frame2 = frame1+framesPerTraj-1;
for jj = 1:numTimeSteps
   if strcmp(dispType,'MFAverage')
      rawDispV{jj} = M(find(M(:,1,frame1)~=0 & M(:,3,frame1)~=0),:,frame1);
      for k = frame1+1:frame2-1
         rawDispV{jj} = [rawDispV{jj}; ...
            M(find(M(:,1,k)~=0 & M(:,3,k)~=0),:,k)];
      end
   elseif strcmp(dispType,'MFTrack')
      MFT          = mFrameTrajBuild(M(:,:,frame1:frame2-1));
      rawDispV{jj} = MFT(:,[1 2 end-1 end]);
   end
   frame1 = frame1+framesPerTStep;
   frame2 = frame1+framesPerTraj-1;

   %Extract raw data points and the raw displamenents on the raw data points.
   rawDataPx{jj} = rawDispV{jj}(:,2);
   rawDataPy{jj} = rawDispV{jj}(:,1);
   rawDataU1{jj} = rawDispV{jj}(:,4)-rawDispV{jj}(:,2);
   rawDataU2{jj} = rawDispV{jj}(:,3)-rawDispV{jj}(:,1);
end

%Plot the raw data points.
plot(rawDataPx{1},rawDataPy{1},'y.'); 

%Load the polygon, 'fieldPG' that define the geometry of the lamollipodium.
load([modelPath 'fieldGeom']);

%Plot the polygon of the lamollipodium.
plot(fieldPGx,fieldPGy,'r-');

%Generate a grid with gridDx X gridDy raster sidelength. And, 
% the grid sampling data points are those inside the predefined polygon.
grid = framework(size(rawI),[gridDy gridDx]);
in   = inpolygon(grid(:,2),grid(:,1),fieldPGx,fieldPGy);
ind  = find(in==1);

gridPx  = cell(numTimeSteps,1);
gridPy  = cell(numTimeSteps,1);
for jj = 1:numTimeSteps
   gridPx{jj} = grid(ind,2);
   gridPy{jj} = grid(ind,1);
end

sDataU1 = cell(numTimeSteps,1);
sDataU2 = cell(numTimeSteps,1);
gridU1  = cell(numTimeSteps,1);
gridU2  = cell(numTimeSteps,1);
if strcmp(calInterp,'yes')
   for jj = 1:numTimeSteps
      %interpolate the vectors on the grid and data points.
      gridUi = vectorFieldInterp(rawDispV{jj},[gridPy{jj} gridPx{jj}], ...
         corLen,[fieldPGx fieldPGy]);
      gridU1{jj} = gridUi(:,4)-gridUi(:,2); %Displacement in x-direction.
      gridU2{jj} = gridUi(:,3)-gridUi(:,1);

      %Smoothed displacements at the raw data points.
      sDispV = vectorFieldInterp(rawDispV{jj},rawDispV(:,1:2), ...
         corLen,[fieldPGx fieldPGy]); 
      sDataU1{jj} = sDispV(:,4)-sDispV(:,2);
      sDataU2{jj} = sDispV(:,3)-sDispV(:,1);
   end
   save([resultPath 'dispField'],'rawI','rawDispV','rawDataPx', ...
      'rawDataPy','rawDataU1','rawDataU2','corLen','edgCorLen', ...
      'gridPx','gridPy','sDataU1','sDataU2','gridU1','gridU2');
else
   save([resultPath 'dispField'],'rawI','rawDispV','rawDataPx', ...
      'rawDataPy','rawDataU1','rawDataU2','corLen','edgCorLen', ...
      'gridPx','gridPy');
end

%Plot the raw data points and displacements.
quiver(rawDataPx{1},rawDataPy{1},rawDataU1{1}*5,rawDataU2{1}*5,0,'y'); 

if strcmp(showInterp,'yes')
   %the filtered displacements.
   quiver(rawDataPx{1},rawDataPy{1},sDataU1{1}*5,sDataU2{1}*5,0,'r'); 
end

hold off;
