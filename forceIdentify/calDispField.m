%Preprocessing of the raw experimental data. It creates by interpolation of
% the raw data a displacement field on grid points inside a predefined polygon.

%Load parameters that are needed for the calculation of the displacement
% field.
run([resultPath 'setPar']);

%Get the image of the cell.
rawI = imread(imgFile);

%Show the image 
figure(gcf); hold off;
imshow(rawI,[]); axis on; hold on;

%Load the displacement field given in the multidimentional array, 'M' that is
% produced by the speckle tracking software.
load(dataFile);

rawDispV = M(find(M(:,1,frameID(1))~=0&M(:,3,frameID(1))~=0),:,frameID(1));
if strcmp(dispType,'MFAverage')
   for k = 2:length(frameID)
      rawDispV = [rawDispV; ...
         M(find(M(:,1,frameID(k))~=0&M(:,3,frameID(k))~=0),:,frameID(k))];
   end
elseif strcmp(dispType,'MFTrack')
   MFT      = mFrameTrajBuild(M(:,:,frameID),length(frameID)+1);
   rawDispV = MFT{1}(:,[1 2 end-1 end]);
end

%Raw data points and the raw displamenents on the raw data points.
rawDataPx = rawDispV(:,2);
rawDataPy = rawDispV(:,1);
rawDataU1 = rawDispV(:,4)-rawDispV(:,2);
rawDataU2 = rawDispV(:,3)-rawDispV(:,1);

%Plot the raw data points.
plot(rawDispV(:,2),rawDispV(:,1),'y.'); 

%Load the polygon, 'fieldPG' that define the geometry of the lamollipodium.
load([modelPath 'fieldGeom']);

%Plot the polygon of the lamollipodium.
plot(fieldPGx,fieldPGy,'r-');

%Generate a grid with gridDx X gridDy raster sidelength. And, 
% the grid sampling data points are those inside the predefined polygon.
grid = framework(size(rawI),[gridDy gridDx]);
in   = inpolygon(grid(:,2),grid(:,1),fieldPGx,fieldPGy);
ind  = find(in==1);

gridPx = grid(ind,2);
gridPy = grid(ind,1);

if strcmp(calInterp,'yes')
   %interpolate the vectors on the grid and data points.
   gridUi = vectorFieldInterp(rawDispV,[gridPy gridPx],corLen,[fieldPGx fieldPGy]);
   gridU1 = gridUi(:,4)-gridUi(:,2); %Displacement in x-direction.
   gridU2 = gridUi(:,3)-gridUi(:,1);

   %Smoothed displacements at the raw data points.
   sDispV  = vectorFieldInterp(rawDispV,rawDispV(:,1:2),corLen,[fieldPGx fieldPGy]); 
   sDataU1 = sDispV(:,4)-sDispV(:,2);
   sDataU2 = sDispV(:,3)-sDispV(:,1);
   save([resultPath 'dispField'],'rawI','rawDispV','rawDataPx', ...
      'rawDataPy','rawDataU1','rawDataU2','corLen','gridPx','gridPy', ...
      'sDataU1','sDataU2','gridU1','gridU2');
else
   save([resultPath 'dispField'],'rawI','rawDispV','rawDataPx', ...
      'rawDataPy','rawDataU1','rawDataU2','corLen','gridPx','gridPy');
end

%Plot the raw data points and displacements.
quiver(rawDataPx,rawDataPy,rawDataU1*5,rawDataU2*5,0,'y'); 

if strcmp(showInterp,'yes')
   %the filtered displacements.
   quiver(rawDataPx,rawDataPy,sDataU1*5,sDataU2*5,0,'r'); 
end

hold off;
