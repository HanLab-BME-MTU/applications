%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Select the polygon that defines the region of interest.
% Always start with the left-lower corner when the polygon is viewed as a
% deformed rectangle and go clockwise.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%choice = input(sprintf(['Draw the field boundary (clockwise) for time step %d:\n' ...
%   '  ''0'': keep the old choice.\n' ...
%   '  ''1'': draw a new polygon.\n' ...
%   'Your answer: '],jj));

if strcmp(isFieldBndFixed,'yes')
   selTimeSteps = 1; %Number of potential field boundaries.
else
   answer = input('Select time steps (0 for all):');
   if isempty(answer) | answer == 0
      selTimeSteps = 1:numDTimePts;
   else
      selTimeSteps = answer;
   end
end

%rawDispFieldDir = [mechDir filesep 'rawDispField'];

figH = [];
numFieldBndsDrawn = 0;
for ii = 1:length(selTimeSteps)
   jj = selTimeSteps(ii);

   imgIndex = imgIndexOfDTimePts(jj);
   rawDispFieldFileName = ['rawDispField' sprintf(imgIndexForm,imgIndex) '.mat'];
   rawDispFieldFile     = [rawDispFieldDir filesep rawDispFieldFileName];
   if exist(rawDispFieldFile,'file')
      s = load(rawDispFieldFile);
      rawDispField = s.rawDispField;
      %rawDataP = s.rawDataP;
      %rawDispV = s.rawDispV;   
   else
      fprintf(1,['Raw vector field has not been loaded and saved for time step %d yet.\n' ...
         'Run loadRawField first.\n'],jj);
      return;
   end

   %Get the image of the cell.
   relFrameNo  = imgIndexOfDTimePts(jj)+relDispImgFrmNo-firstImgIndex;
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

   %Show the raw displacement.
   quiver(rawDispField.p(rawDispField.inlier,1),rawDispField.p(rawDispField.inlier,2), ...
      rawDispField.v(rawDispField.inlier,1)*dispScale,rawDispField.v(rawDispField.inlier,2)*dispScale,0,'y');
   quiver(rawDispField.p(rawDispField.outlier,1),rawDispField.p(rawDispField.outlier,2), ...
      rawDispField.v(rawDispField.outlier,1)*dispScale,rawDispField.v(rawDispField.outlier,2)*dispScale,0,'r');
   %quiver(rawDispV{1}(rawInlier{1},2),rawDispV{1}(rawInlier{1},1), ...
   %   (rawDispV{1}(rawInlier{1},4)-rawDispV{1}(rawInlier{1},2))*dispScale, ...
   %   (rawDispV{1}(rawInlier{1},3)-rawDispV{1}(rawInlier{1},1))*dispScale,0,'y');
   %quiver(rawDispV{1}(rawOutlier{1},2),rawDispV{1}(rawOutlier{1},1), ...
   %   (rawDispV{1}(rawOutlier{1},4)-rawDispV{1}(rawOutlier{1},2))*dispScale, ...
   %   (rawDispV{1}(rawOutlier{1},3)-rawDispV{1}(rawOutlier{1},1))*dispScale,0,'r');

   if numFieldBndsDrawn >= 1
      plot(fieldPGx,fieldPGy,'w-.');
   end

   %Identify the file where the boundary is stored.
   if strcmp(isFieldBndFixed,'yes')
      titleStr = sprintf(['Draw the field boundary for all time steps\n' ...
         'White: previous drawn; Green: last saved']);
      qstr = '';
      femModelFile = [femModelDir filesep 'femModel' sprintf(imgIndexForm,0) '.mat'];
      numFieldBndsDrawn = 1;
   else
      qstr = sprintf('Do you want to draw the field boundary for time step: %d?(y/n)',jj); 
      titleStr = sprintf(['Draw the field boundary for time step: %d; Image Index: %d\n' ...
         'White: previous time step; Green: last saved'],jj,imgIndex);
      femModelFile = [femModelDir filesep 'femModel' sprintf(imgIndexForm,imgIndex) '.mat'];
   end

   if exist(femModelFile,'file') == 2
      s = load(femModelFile);
      femModel = s.femModel;

      plot(femModel.fieldBnd.x,femModel.fieldBnd.y,'g-.');
   end

   title(titleStr);

   if ~isempty(qstr)
      ans = input(qstr,'s');
   else
      ans = 'y';
   end

   if strcmp(ans,'y')

      numFieldBndsDrawn = numFieldBndsDrawn+1;
      %Draw the boundary of the field.
      [bw,fieldPGx,fieldPGy] = roipoly;
      fieldBnd.x = fieldPGx;
      fieldBnd.y = fieldPGy;

      plot(fieldPGx,fieldPGy,'r');
      for k = 1:length(fieldPGx)
         text(fieldPGx(k),fieldPGy(k),num2str(k),'color','g');
      end

      %Ask for the indices of the other three vertices of the polygon field boundary 
      % when it is viewed as a deformed rectangle. The first vertex is always 1.
      fieldPGVI = ones(1,4);
      fieldPGVI(2:4) = input(['What are the indices of the four vertices ' ...
         'when the polygon is viewed as a deformed rectangle?\n  ' ...
         'Start with the left-lower corner (index: 1) and\n  ' ...
         'enter the other three by going clockwise:']);
      plot(fieldPGx(fieldPGVI),fieldPGy(fieldPGVI),'go');

      fieldBnd.VI = fieldPGVI;
      femModel.fieldBnd = fieldBnd;

      plot(femModel.fieldBnd.x,femModel.fieldBnd.y,'r-.');
      save(femModelFile,'femModel');
   end
   %save([mechDir filesep 'fieldGeom'], 'fieldPGx', 'fieldPGy');
end

%for k = 1:length(fieldPGx)
%   text(fieldPGx(k),fieldPGy(k),num2str(k),'color','r');
%end
%plot(fieldPGx,fieldPGy,'go');
%plot(fieldPGx,fieldPGy,'b');
