%Preprocessing of the raw experimental data. It creates by interpolation of
% the raw data a displacement field for force reconstruction.

%rawDispFieldDir = [mechDir filesep 'rawDispField'];

%Generate a grid with gridDx X gridDy raster sidelength. 
[numPixelsY,numPixelsX] = size(firstDTImg);
numGridX = ceil(numPixelsX/gridDx);
numGridY = ceil(numPixelsY/gridDy);
gridx    = linspace(1,numPixelsX,numGridX);
gridy    = linspace(1,numPixelsY,numGridY);
[gridX,gridY] = meshgrid(gridx,gridy);
gridX = reshape(gridX,length(gridX(:)),1);
gridY = reshape(gridY,length(gridY(:)),1);

figH = [];
for jj = 1:numDTimePts
   imgIndex = imgIndexOfDTimePts(jj);
   %Identify the file where the raw displacement field is stored.
   rawDispFieldFileName = ['rawDispField' sprintf(imgIndexForm,imgIndex) '.mat'];
   rawDispFieldFile     = [rawDispFieldDir filesep rawDispFieldFileName];

   %Load the raw field
   if exist(rawDispFieldFile,'file')
      s = load(rawDispFieldFile);
      rawDispField = s.rawDispField;
   else
      fprintf(1,'Raw vector field has not been loaded and saved yet. Run loadRawField first.\n');
      return;
   end

   %Identify the file where the boundary is stored.
   if strcmp(isFieldBndFixed,'yes')
      femModelFile = [femModelDir filesep 'femModel' sprintf(imgIndexForm,0) '.mat'];
   else
      femModelFile = getFemModelFile(femModelDir,imgIndex,imgIndexForm);
   end

   %Load the boundary.
   if exist(femModelFile,'file')
      s = load(femModelFile);
      fieldBnd = s.femModel.fieldBnd;
   else
      fprintf(1,'Field boundary has not been drawn and saved yet. Run drawFieldBound first.\n');
      return;
   end

   %Get the image of the cell.
   relFrameNo = imgIndexOfDTimePts(jj)+relDispImgFrmNo-firstImgIndex;
   overlaidImg = double(imread(imgFileList{1}{relFrameNo}));
   for kk = 1:numAvgFrames(jj)-1
      overlaidImg = overlaidImg+double(imread(imgFileList{1}{relFrameNo+kk}));
   end
   overlaidImg = overlaidImg/numAvgFrames(jj);

   %Show the image 
   if isempty(figH) || ~ishandle(figH)
      figH = figure;
   end
   figure(figH); hold off;
   imshow(overlaidImg,[]); axis on; hold on;

   title(sprintf('Calculating the displacement field for time step: %d',jj));

   %Identify grid points that are inside the field boundary.
   in   = inpolygon(gridX,gridY,fieldBnd.x,fieldBnd.y);
   ind  = find(in==1);

   gridPx = gridX(ind);
   gridPy = gridY(ind);


   %We interpolate raw displacements to [iDataPx iDataPy] with very small
   % space correlation length: 'MFCorLen'.
   %iDataPx = cell(numDTimePts,1);
   %iDataPy = cell(numDTimePts,1);
   %iDataU1 = cell(numDTimePts,1);
   %iDataU2 = cell(numDTimePts,1);

   %We interpolate raw displacements to points which can be either the original data
   % position or the grid points.
   rawDispV = [rawDispField.p(:,2:-1:1) rawDispField.p(:,2:-1:1)+rawDispField.v(:,2:-1:1)];
   if strcmp(trackMethod,'speck')
      frame1 = imgIndexOfDTimePts(jj)-firstImgIndex+1;
      frame2 = frame1+numAvgFrames-1;
      if strcmp(dataSite,'grid')
         iDispField.p = [gridPx gridPy];
      elseif strcmp(dataSite,'original')
         in  = inpolygon(rawDispField.p(:,1),rawDispField.p(:,2),fieldBnd.x,fieldBnd.y);
         ind = find(in==1);

         iDispField.p = rawDispField.p(ind,:);
      elseif strcmp(dataSite,'everyOther')
         in  = inpolygon(rawDispField.p(:,1),rawDispField.p(:,2),fieldBnd.x,fieldBnd.y);
         ind = find(in==1);

         iDispField.p = rawDispField.p(ind(1:2:end),:);
      else
         error(['''dataSite'': ' dataSite ' is not recogonized.']);
      end

      if strcmp(dispType,'MFAverage')
         dataDispV = vectorFieldSparseInterp(rawDispV, ...
            iDispField.p(:,2:-1:1),2*corLen,corLen,[]); 
      elseif strcmp(dispType,'MFTrack')
         %MFT          = mFrameTrajBuild(M(:,:,frame1:frame2-1));
         MFT = mFrameTrajBuild(M(:,:,frame1:frame2-1),corLen, ...
            iDispField.p(:,2:-1:1));
         dataDispV = MFT(:,[1 2 end-1 end]);
      elseif strcmp(dispType,'SFrame')
         if strcmp(dataSite,'grid')
            dataDispV = vectorFieldSparseInterp(rawDispV, ...
               iDispField.p(:,2:-1:1),2*corLen,corLen,[]); 
         else
            dataDispV = rawDispV;
         end
      end

      %Extract data points and the raw displamenents on the data points.
      iDispField.p = dataDispV(:,2:-1:1);
      iDispField.v = dataDispV(:,4:-1:3)-dataDispV(:,2:-1:1);
   elseif strcmp(trackMethod,'corr') || strcmp(trackMethod,'dArray')
      if strcmp(dataSite,'grid')
         iDispField.p = [gridPx gridPy];
         dataDispV = vectorFieldSparseInterp(rawDispV, ...
            iDispField.p(:,2:-1:1),2*corLen,corLen,[]); 
         iDispField.v = dataDispV(:,4:-1:3)-dataDispV(:,2:-1:1);
      elseif strcmp(dataSite,'original')
         in  = inpolygon(rawDispField.p(:,1),rawDispField.p(:,2),fieldBnd.x,fieldBnd.y);
         ind = find(in==1);

         iDispField.p = rawDispField.p(ind,:);
         iDispField.v = rawDispField.v(ind,:);
      elseif strcmp(dataSite,'everyOther')
         in  = inpolygon(rawDispField.p(:,1),rawDispField.p(:,2),fieldBnd.x,fieldBnd.y);
         ind = find(in==1);

         iDispField.p = rawDispField.p(ind(1:2:end),:);
         iDispField.v = rawDispField.v(ind(1:2:end),:);
      end
   else
      error(['''trackMethod'': ' trackMethod ' is not recogonized.']);
   end

   nanInd = find(isnan(iDispField.v(:,1)) | isnan(iDispField.v(:,2)));
   iDispField.p(nanInd,:) = [];
   iDispField.v(nanInd,:) = [];

   %Plot the data points and interpolated displacements.
   quiver(iDispField.p(:,1),iDispField.p(:,2),iDispField.v(:,1)*dispScale, ...
      iDispField.v(:,2)*dispScale,0,'y'); 

   %We also calculate the interpolated displacements with a bigger correlation length.
   % They are considered smoothed displacements at the data points.
   sDispV = vectorFieldSparseInterp(rawDispV, ...
      iDispField.p(:,2:-1:1), 2*sCorLen,sCorLen,[fieldBnd.x fieldBnd.y]); 
   iDispField.sv = sDispV(:,4:-1:3)-sDispV(:,2:-1:1);

   if strcmp(showSmooth,'yes')
      %the filtered displacements.
      quiver(iDispField.p(:,1),iDispField.p(:,2),iDispField.sv(:,1)*dispScale, ...
         iDispField.sv(:,2)*dispScale,0,'r'); 
   end

   iDispField.corLen     = corLen;
   iDispField.sCorLen    = sCorLen;
   iDispField.edgCorLen  = edgCorLen;
   iDispField.gridx      = gridx;
   iDispField.gridy      = gridy;
   iDispField.gridX      = gridX;
   iDispField.gridY      = gridY;
   iDispField.gridPx     = gridPx;
   iDispField.gridPy     = gridPy;
   %Identify the file where the interpolated displacement field is stored.
   iDispFieldFileName = ['iDispField' sprintf(imgIndexForm,imgIndex) '.mat'];
   iDispFieldFile     = [iDispFieldDir filesep iDispFieldFileName];

   save(iDispFieldFile,'iDispField');
   %save([reslDir filesep 'iDispField'],'iDataPx', ...
   %   'iDataPy','iDataU1','iDataU2','sDataU1','sDataU2', ...
   %   'MFCorLen','corLen','edgCorLen','gridPx','gridPy', ...
   %   'gridx','gridy','gridX','gridY');
end

