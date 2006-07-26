%Post assemble of the identified body force and some calculation
% for debugging purpose. It is run after 'solveLSBF'.

startTime = cputime;

fprintf(1,'Post assemble of cellular domain force:\n');

%Get some physical parameters.
load([projDir filesep 'lastProjSettings.mat']);
physiParam = projSettings.physiParam;

actPixelSize  = 67; %Unit, nm.
frameInterval = 10; %Unit, sec.
if iscell(physiParam)
   %We assume the actin is the 1st image channel.
   actPixelSize  = physiParam{1}.pixelSize;
   frameInterval = physiParam{1}.frameInterval;
end

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

answer = input('Select time steps (0 for all):');
if isempty(answer) | answer == 0
   selTimeSteps = 1:numDTimePts;
else
   selTimeSteps = answer;
end

for ii = 1:length(selTimeSteps)
   jj = selTimeSteps(ii);

   localStartTime = cputime;

   fprintf(1,'   Time step %d ... ',jj);

   imgIndex = imgIndexOfDTimePts(jj);

   %Get one cell image
   dispImgIndex  = max(firstImgIndex,imgIndex) + relDispImgFrmNo - 1;
   curRelFrmNo   = dispImgIndex-firstImgIndex+1;
   dispImg = imread(imgFileList{imgChannel}{curRelFrmNo});

   %Get the stacked cell image in the averaging range.
   frameNo = max(firstImgIndex,imgIndex)-firstImgIndex+1;
   if length(numAvgFrames) == numDTimePts
      curNumAvgFrames = numAvgFrames(curDTimePt);
   else
      curNumAvgFrames = numAvgFrames;
   end

   %Get the overlaid images of selected image channel.
   stackedImg = double(imread(imgFileList{imgChannel}{frameNo}));
   for k = 1:curNumAvgFrames-1
      stackedImg = stackedImg + double(imread(imgFileList{imgChannel}{k+frameNo}));
   end
   stackedImg = stackedImg/curNumAvgFrames;

   %if strcmp(bfFwdOpComputed,'all') == 1
   %   fprintf(1,'Load the identified body force and displacements: \n');
   %   load([reslDir filesep 'bfId']);
   %   load([reslDir filesep 'dispId']);
   %   return;
   %end

   if strcmp(isFieldBndFixed,'yes')
      femModelFile = [femModelDir filesep 'femModel' sprintf(imgIndexForm,0) '.mat'];
   else
      [femModelFile,femModelImgIndex] = getFemModelFile(femModelDir,imgIndex,imgIndexForm);
   end
   s = load(femModelFile);
   femModel = s.femModel;
   fem      = femModel.fem;
   numEdges = femModel.numEdges;
   fn       = femModel.fn;
   fp       = femModel.fp;
   fs       = femModel.fs;

   %Specify boundary condition.
   for k = 1:numEdges
      fn.BndDispx{k} = 'bndDisp';
      fn.BndDispy{k} = 'bndDisp';
   end

   fn.BodyFx = 'femBodyF';
   fn.BodyFy = 'femBodyF';

   iDispFieldFileName = ['iDispField' sprintf(imgIndexForm,imgIndex) '.mat'];
   iDispFieldFile     = [iDispFieldDir filesep iDispFieldFileName];
   s = load(iDispFieldFile);
   iDispField = s.iDispField;

   rawDispFieldFileName = ['rawDispField' sprintf(imgIndexForm,imgIndex) '.mat'];
   rawDispFieldFile     = [rawDispFieldDir filesep rawDispFieldFileName];
   s = load(rawDispFieldFile);
   rawDispField = s.rawDispField;

   %For debugging : Solve the elastic equation with the identified force and
   % compare the computed displacement with the measured displacement data.

   %Get the mask of cell from 'edgeDir'. It will be used in creating data map.
   cellMask = [];
   if isdir(edgeDir)
      maskDir = [edgeDir filesep 'cell_mask'];
      if isdir(maskDir)
         [maskFileList maskIndex] = getNamedFiles(maskDir,'mask_');

         thisMaskIndex = find(maskIndex==dispImgIndex);
         if ~isempty(thisMaskIndex)
            cellMaskFile = maskFileList{thisMaskIndex(1)};
            cellMask = imread([maskDir filesep cellMaskFile]);
         end
      end
   end

   dimBF     = fs.dimBF;
   dimFS     = fs.dimFS;
   indDomDOF = fs.indDomDOF;

   forceFieldFile = [forceFieldDir filesep 'forceField' ...
      sprintf(imgIndexForm,imgIndex) '.mat'];
   s = load(forceFieldFile);
   forceField = s.forceField;

   coefBFx = forceField.coefBF(:,1);
   coefBFy = forceField.coefBF(:,2);

   fp.BodyFx = {{'x' 'y'} {fs.fem coefBFx}};
   fp.BodyFy = {{'x' 'y'} {fs.fem coefBFy}};

   edgD = iDispField.edgD;
   for k = 1:numEdges
      fp.BndDispx{k} = {{'s'} {edgD(k).ppU1}};
      fp.BndDispy{k} = {{'s'} {edgD(k).ppU2}};
   end
   %for k = 1:numEdges
   %   fp.BndDispx{k} = {{'s'} {edgePPx{jj,k}}};
   %   fp.BndDispy{k} = {{'s'} {edgePPy{jj,k}}};
   %end

   fem = elModelUpdate(fem,'fn',fn,'fp',fp);
   fem = elasticSolve(fem,[]);

   %The data points must also be inside the target recovery region.
   forceField.p = iDispField.p(iDispField.iInMesh,:);
   %[is,pe] = postinterp(fem,forceField.p.');
   %forceField.p(pe,:) = [];

   %Calculate the identified force.
   %[recBFx,recBFy] = postinterp(fem,'f1','f2', ...
   %   [bfDisplayPx bfDisplayPy].');
   %forceField.recBF = [recBFx;recBFy].';
   %recBF(:,:,jj) = [recBFx; recBFy].';
   %[bfxR,bfyR] = postinterp(fem,'f1','f2',fs.mesh.p);
   [recBFx,recBFy] = postinterp(fem,'f1','f2', ...
      forceField.p.');
   forceField.f = [recBFx;recBFy].';

   %The displacements on the points where the force field is to be
   % demonstrated.
   %[recDispU1,recDispU2] = postinterp(fem,'u1','u2', ...
   %   [bfDisplayPx bfDisplayPy].');
   [recDispU1,recDispU2] = postinterp(fem,'u1','u2', ...
      forceField.p.');
   iDispField.rv = [recDispU1;recDispU2].';
   %recDispU(:,:,jj) = [recDispU1;recDispU2].';

   %Get the grid points that are inside the identification region. They will be
   % used to calculate the intensity score for contraction and adhesion.
   gridx = iDispField.gridx;
   gridy = iDispField.gridy;
   gridX = iDispField.gridX;
   gridY = iDispField.gridY;
   [is,pe] = postinterp(fem,[gridX gridY].');
   gridIn     = 1:length(gridX);
   gridIn(pe) = [];
   gridXin    = gridX(gridIn);
   gridYin    = gridY(gridIn);

   forceField.gridIn = gridIn;

   %Calculate the strain.
   %[recStrnU11,recStrnU12,recStrnU22] = postinterp(fem,'u1x', ...
   %   '(u1y+u2x)/2','u2y',[gridX(gridIn) gridY(gridIn)].');
   %recStrain(gridIn,:,jj) = [recStrnU11;recStrnU12;recStrnU22].';
   [recStrnU11,recStrnU12,recStrnU22] = postinterp(fem,'u1x', ...
      '(u1y+u2x)/2','u2y',[gridXin gridYin].');
   recStrain = [recStrnU11;recStrnU12;recStrnU22].';

   %Calculate the stress.
   %[dispLambda,dispMu] = postinterp(fem,'lambda','mu', ...
   %   [gridX(gridIn) gridY(gridIn)].');
   [gridLambda,gridMu] = postinterp(fem,'lambda','mu', ...
      [gridXin gridYin].');
   recStress = (gridLambda.*(recStrnU11+recStrnU22) + ...
      2*gridMu.*recStrnU11).';
   recStress = [recStress; (2*gridMu.*recStrnU12).'];
   recStress = [recStress; (gridLambda.*(recStrnU11+recStrnU22) + ...
      2*gridMu.*recStrnU22).'];

   %for k = 1:3
   %   ppStress{k,jj} = csape({gridy,gridx}, ...
   %      reshape(squeeze(recStress(:,k,jj)),length(gridy),length(gridx)));
   %   ppStrain{k,jj} = csape({gridy,gridx}, ...
   %      reshape(squeeze(recStrain(:,k,jj)),length(gridy),length(gridx)));
   %end

   %%%%%%% Identify the location of adhesion and contration and assign intensity
   % scores that approximately separate the two. It is based on the angle
   % between the force and the displacement.
   %We first separate the two for forces on the display points.
   recBFLen = sqrt(recBFx.^2+recBFy.^2);
   recDispLen = sqrt(recDispU1.^2+recDispU2.^2);

   unitRecU1 = recDispU1./recDispLen;
   unitRecU2 = recDispU2./recDispLen;

   dotProdBFRecU = recBFx.*unitRecU1 + recBFy.*unitRecU2;
   mcfIndShow = find((dotProdBFRecU-recBFLen*cos(mcfAngle))>=0);
   adfIndShow = find((-dotProdBFRecU-recBFLen*cos(adfAngle))>=0);
   mixIndShow = find((dotProdBFRecU-recBFLen*cos(mcfAngle))<0 & ...
       (-dotProdBFRecU-recBFLen*cos(adfAngle))<0);

   recADFx = NaN*ones(size(recBFx));
   recADFy = recADFx;
   recADFx(adfIndShow) = recBFx(adfIndShow);
   recADFy(adfIndShow) = recBFy(adfIndShow);
    
   recMCFx = NaN*ones(size(recBFx));
   recMCFy = recMCFx;
   recMCFx(mcfIndShow) = recBFx(mcfIndShow);
   recMCFy(mcfIndShow) = recBFy(mcfIndShow);
   
   %%%%%%% For mixed zone, we project the force into the opposite direction of
   % flow (adhesion force) and the direction in between the force vector and 
   % the vector that forms 'mcfAngle' with the flow (contraction force).
   % First, get the direction vector that forms 'mcfAngle' with the flow.
   % This is nothing but a rotation of the flow vector towards the force
   % vector by 'mcfAngle'. The sign of rotation angle is given by the sign
   % of the dot product between 'force' and the orthogonal vector to the
   % flow.
   dotProdBF_RecUPerp = recBFx.*unitRecU2 - recBFy.*unitRecU1;
   
   %'(mcfAnglV1,mcfAnglV2)': the unit vector that forms 'mcfAngle' with
   % the unit flow vector.
   rotAngle = sign(dotProdBF_RecUPerp(mixIndShow))*mcfAngle;
   mcfAnglV1 = cos(rotAngle).*unitRecU1(mixIndShow) + ...
       sin(rotAngle).*unitRecU2(mixIndShow);
   mcfAnglV2 = -sin(rotAngle).*unitRecU1(mixIndShow) + ...
       cos(rotAngle).*unitRecU2(mixIndShow);
   
   if strcmp(fixMcfProj,'yes')
      mcfUnitV1 = mcfAnglV1;
      mcfUnitV2 = mcfAnglV2;
   else
      %'(mcfUnitV1,mcfUnitV2)': the unit vector in between '(mcfAnglV1,mcfAnglV2)'
      % and the force vector.
      mcfUnitV1   = recBFx(mixIndShow)./recBFLen(mixIndShow) + mcfAnglV1;
      mcfUnitV2   = recBFy(mixIndShow)./recBFLen(mixIndShow) + mcfAnglV2;
      mcfUnitVLen = sqrt(mcfUnitV1.^2 + mcfUnitV2.^2);
      mcfUnitV1   = mcfUnitV1./mcfUnitVLen;
      mcfUnitV2   = mcfUnitV2./mcfUnitVLen;
   end

   %'(adfAnglV1,adfAnglV2)': the unit vector that forms 'adfAngle' with
   % the unit flow vector.
   rotAngle = -sign(dotProdBF_RecUPerp(mixIndShow))*adfAngle;
   adfUnitV1 = -cos(rotAngle).*unitRecU1(mixIndShow) - ...
       sin(rotAngle).*unitRecU2(mixIndShow);
   adfUnitV2 = sin(rotAngle).*unitRecU1(mixIndShow) - ...
       cos(rotAngle).*unitRecU2(mixIndShow);
   
   for k = 1:length(mixIndShow)
       %'projM': The two columns of 'projM' are given by the vector
       %opposite the flow and the unit vector ((mcfUnitV1,mcfUnitV2) that
       %forms 'mcfAngle' with the flow.
       %projM = [-unitRecU1(mixIndShow(k)) mcfUnitV1(k); ...
       %        -unitRecU2(mixIndShow(k)) mcfUnitV2(k)];
       projM = [adfUnitV1(k) mcfUnitV1(k); ...
               adfUnitV2(k) mcfUnitV2(k)];
       %The amount of projection to the two vectors are given by
       pC = projM\[recBFx(mixIndShow(k));recBFy(mixIndShow(k))];
       recADFx(mixIndShow(k)) = projM(1,1)*pC(1);
       recADFy(mixIndShow(k)) = projM(2,1)*pC(1);
       recMCFx(mixIndShow(k)) = projM(1,2)*pC(2);
       recMCFy(mixIndShow(k)) = projM(2,2)*pC(2);
   end
   forceField.adf = [recADFx; recADFy].';
   forceField.mcf = [recMCFx; recMCFy].';
   
   %%%%%%% We then do it on the grid points and assign intensity scores.
   % The force on grid points inside the identification region.
   %We calculate forces on grids for generating color map.
   gridBDF  = NaN*ones(length(gridy)*length(gridx),1);
   gridMCF  = gridBDF;
   gridADF  = gridBDF;
   gridADC  = gridBDF;
   gridMixF = gridBDF;
   gridSpd  = gridBDF;

   %We can also display the used Young's modulus.
   gridYMod = gridBDF;

   [gridBFx,gridBFy] = postinterp(fem,'f1','f2',[gridXin gridYin].');
   forceField.gridF = [gridBFx; gridBFy].';

   gridInMCFx  = NaN*ones(size(gridBFx));
   gridInMCFy  = NaN*ones(size(gridBFy));
   gridInADFx  = NaN*ones(size(gridBFx));
   gridInADFy  = NaN*ones(size(gridBFy));
   gridInMixFx = NaN*ones(size(gridBFx));
   gridInMixFy = NaN*ones(size(gridBFy));

   %The displacements on grid points.
   [gridU1,gridU2] = postinterp(fem,'u1','u2',[gridXin gridYin].');
   forceField.gridV = [gridU1;gridU2].';

   %Young's modulus on grid points inside the identification region.
   gridInYMod = postinterp(fem,'YModul',[gridXin gridYin].');
   gridYMod(gridIn) = gridInYMod;

   %Compute the length of the body force.
   gridBFLen = sqrt(gridBFx.^2+gridBFy.^2);

   %Compute the length of the displacements and their average .
   gridDispLen = sqrt(gridU1.^2+gridU2.^2);
   avgDispLen  = sum(gridDispLen(:))/length(gridDispLen(:));
   smDispLen   = avgDispLen*smDispThreshold;

   %Normalize the displacement.
   %unitGridU1 = zeros(size(gridU1));
   %unitGridU2 = zeros(size(gridU2));
   %unitGridU1(sigDispInd) = gridU1(sigDispInd)./gridDispLen(sigDispInd);
   %unitGridU2(sigDispInd) = gridU2(sigDispInd)./gridDispLen(sigDispInd);
   unitGridU1 = gridU1./gridDispLen;
   unitGridU2 = gridU2./gridDispLen;
   
   %Calculate the intensity scores for contraction and adhesion.
   % 'gridInMCF' : Score for Myosin Contractile Force.
   % 'gridInADF' : Score for Adhesion Dragging Force.
   %The sign of the dot product between the force and the displacement tells
   % if the angle between the force and the displacment is less or greater than
   % pi/2.
   %gridInADF = -gridBFLen;
   %gridInADF(sigDispInd) = gridBFx(sigDispInd).*unitGridU1(sigDispInd) + ...
   %   gridBFy(sigDispInd).*unitGridU2(sigDispInd);
   %First comput the dot product between the force and the displacements.
   dotProdBFU = gridBFx.*unitGridU1 + gridBFy.*unitGridU2;

   % mcfInd : Indices where we consider the total body force to be mainly 
   %          myosin contraction.
   % adfInd : Indices where we consider the total body force to be mainly 
   %          adhesion resistance.
   % mixInd : Indices where we consider the total body force to be a mixture 
   %          of myosin contraction and adhesion resistance.
   mcfInd = find((dotProdBFU-gridBFLen*cos(mcfAngle))>=0);
   adfInd = find((-dotProdBFU-gridBFLen*cos(adfAngle)>=0));
   mixInd = find(dotProdBFU-gridBFLen*cos(mcfAngle)<0 & ...
      -dotProdBFU-gridBFLen*cos(adfAngle)<0);

   gridInADF  = NaN*ones(size(gridBFLen));
   gridInADC  = NaN*ones(size(gridBFLen)); %Adhesion drag coefficient.
   gridInMCF  = gridInADF;
   gridInMixF = gridInADF;

   %Dominant myosin dragging force.
   gridInMCF(mcfInd)  = gridBFLen(mcfInd);
   gridInMCFx(mcfInd) = gridBFx(mcfInd);
   gridInMCFy(mcfInd) = gridBFy(mcfInd);
   
   %Dominant adhesion force.
   gridInADF(adfInd)  = gridBFLen(adfInd);
   gridInADFx(adfInd) = gridBFx(adfInd);
   gridInADFy(adfInd) = gridBFy(adfInd);

   %Mixture of myosin and adhesion force.
   gridInMixF(mixInd)  = gridBFLen(mixInd);
   gridInMixFx(mixInd) = gridBFx(mixInd);
   gridInMixFy(mixInd) = gridBFy(mixInd);

   %%%%%%% For mixed zone, we project the force into the opposite direction of
   % flow (adhesion force) and the direction that forms 'mcfAngle' with 
   % the flow (contraction force).
   % First, get the direction vector that forms 'mcfAngle' with the flow.
   % This is nothing but a rotation of the flow vector towards the force
   % vector by 'mcfAngle'. The sign of rotation angle is given by the sign
   % of the dot product between 'force' and the orthogonal vector to the
   % flow.
   dotProdBF_GridUPerp = gridBFx.*unitGridU2 - gridBFy.*unitGridU1;
   
   %'(mcfAnglV1,mcfAnglV2)': the unit vector that forms 'mcfAngle' with
   % the unit flow vector.
   rotAngle = sign(dotProdBF_GridUPerp(mixInd))*mcfAngle;
   mcfAnglV1 = cos(rotAngle).*unitGridU1(mixInd) + ...
       sin(rotAngle).*unitGridU2(mixInd);
   mcfAnglV2 = -sin(rotAngle).*unitGridU1(mixInd) + ...
       cos(rotAngle).*unitGridU2(mixInd);
   
   if strcmp(fixMcfProj,'yes')
      mcfUnitV1 = mcfAnglV1;
      mcfUnitV2 = mcfAnglV2;
   else
      %'(mcfUnitV1,mcfUnitV2)': the unit vector in between '(mcfAnglV1,mcfAnglV2)'
      % and the force vector.
      mcfUnitV1   = gridBFx(mixInd)./gridBFLen(mixInd) + mcfAnglV1;
      mcfUnitV2   = gridBFy(mixInd)./gridBFLen(mixInd) + mcfAnglV2;
      mcfUnitVLen = sqrt(mcfUnitV1.^2 + mcfUnitV2.^2);
      mcfUnitV1   = mcfUnitV1./mcfUnitVLen;
      mcfUnitV2   = mcfUnitV2./mcfUnitVLen;
   end

   %'(adfAnglV1,adfAnglV2)': the unit vector that forms 'adfAngle' with
   % the unit flow vector.
   rotAngle = -sign(dotProdBF_GridUPerp(mixInd))*adfAngle;
   adfUnitV1 = -cos(rotAngle).*unitGridU1(mixInd) - ...
       sin(rotAngle).*unitGridU2(mixInd);
   adfUnitV2 = sin(rotAngle).*unitGridU1(mixInd) - ...
       cos(rotAngle).*unitGridU2(mixInd);

   for k = 1:length(mixInd)
       %'projM': The two columns of 'projM' are given by the vector
       %opposite the flow and the unit vector ((mcfUnitV1,mcfUnitV2) that
       %forms 'mcfAngle' with the flow.
       %projM = [-unitGridU1(mixInd(k)) mcfUnitV1(k); ...
       %        -unitGridU2(mixInd(k)) mcfUnitV2(k)];
       projM = [adfUnitV1(k) mcfUnitV1(k); ...
               adfUnitV2(k) mcfUnitV2(k)];
       %The amount of projection to the two vectors are given by
       pC = projM\[gridBFx(mixInd(k));gridBFy(mixInd(k))];
       gridInADF(mixInd(k))  = pC(1);
       gridInMCF(mixInd(k))  = pC(2);
       gridInADFx(mixInd(k)) = projM(1,1)*pC(1);
       gridInADFy(mixInd(k)) = projM(2,1)*pC(1);
       gridInMCFx(mixInd(k)) = projM(1,2)*pC(2);
       gridInMCFy(mixInd(k)) = projM(2,2)*pC(2);
   end
   gridInADF(mixInd(find(gridInADF(mixInd)<gridInMCF(mixInd)))) = NaN;

   %gridInMCF(mixInd)  = sqrt(gridInMCFx(mixInd).^2+gridInMCFy(mixInd).^2);
  %gridInMCF = zeros(size(gridBFLen));
   %gridInMCF(sigDispInd) = abs(gridBFx(sigDispInd).*unitGridU2(sigDispInd) - ...
   %   gridBFy(sigDispInd).*unitGridU1(sigDispInd));
   %gridInMCF(mcfInd) = gridBFLen(mcfInd);

   %Calculate the drag coefficient. For small displacements that are 
   % below a threshold value, we divide by the threshold.
   % 'smDispInd' : Index of the small displacements.
   smDispInd  = find(gridDispLen<smDispLen);
   bigDispInd = find(gridDispLen>=smDispLen);

   gridInADC(smDispInd)  = gridInADF(smDispInd)./smDispLen;
   gridInADC(bigDispInd) = gridInADF(bigDispInd)./gridDispLen(bigDispInd);

   %Interpolate 'gridInMCF' and 'gridInADF' to get the pp-form of the scores.
   gridBDF(gridIn)  = gridBFLen;
   gridMCF(gridIn)  = gridInMCF;
   gridADF(gridIn)  = gridInADF;
   gridADC(gridIn)  = gridInADC;
   gridMixF(gridIn) = gridInMixF;
   gridSpd(gridIn)  = sqrt(gridU1.^2+gridU2.^2);

   gridBDF  = reshape(gridBDF,length(gridy),length(gridx));
   gridMCF  = reshape(gridMCF,length(gridy),length(gridx));
   gridADF  = reshape(gridADF,length(gridy),length(gridx));
   gridADC  = reshape(gridADC,length(gridy),length(gridx));
   gridMixF = reshape(gridMixF,length(gridy),length(gridx));
   gridYMod = reshape(gridYMod,length(gridy),length(gridx));
   gridSpd  = reshape(gridSpd,length(gridy),length(gridx));

   forceField.gridBDF  = gridBDF;
   forceField.gridMCF  = gridMCF;
   forceField.gridADF  = gridADF;
   forceField.gridADC  = gridADC;
   forceField.gridMixF = gridMixF;
   forceField.gridYMod = gridYMod;
   forceField.gridSpd  = gridSpd;

   if isempty(cellMask)
      bdfMap = imDataMap(size(stackedImg),{gridy,gridx},gridBDF,'bnd',[bfDomPGx bfDomPGy]);
   else
      bdfMap = imDataMap(size(stackedImg),{gridy,gridx},gridBDF,'mask',cellMask);
   end

   if isempty(cellMask)
      mcfMap  = imDataMap(size(stackedImg),{gridy,gridx},gridMCF,'bnd',[bfDomPGx bfDomPGy]);
   else
      mcfMap  = imDataMap(size(stackedImg),{gridy,gridx},gridMCF,'mask',cellMask);
   end

   if isempty(cellMask)
      adfMap  = imDataMap(size(stackedImg),{gridy,gridx},gridADF,'bnd',[bfDomPGx bfDomPGy]);
   else
      adfMap  = imDataMap(size(stackedImg),{gridy,gridx},gridADF,'mask',cellMask);
   end

   if isempty(cellMask)
      adcMap  = imDataMap(size(stackedImg),{gridy,gridx},gridADC,'bnd',[bfDomPGx bfDomPGy]);
   else
      adcMap  = imDataMap(size(stackedImg),{gridy,gridx},gridADC,'mask',cellMask);
   end

   if isempty(cellMask)
      mixfMap = imDataMap(size(stackedImg),{gridy,gridx},gridMixF,'bnd',[bfDomPGx bfDomPGy]);
   else
      mixfMap = imDataMap(size(stackedImg),{gridy,gridx},gridMixF,'mask',cellMask);
   end

   if isempty(cellMask)
      spdMap  = imDataMap(size(stackedImg),{gridy,gridx},gridSpd,'bnd',[bfDomPGx bfDomPGy]);
   else
      spdMap  = imDataMap(size(stackedImg),{gridy,gridx},gridSpd,'mask',cellMask);
   end

   if isempty(cellMask)
      ymodMap = imDataMap(size(stackedImg),{gridy,gridx},gridYMod,'bnd',[bfDomPGx bfDomPGy]);
   else
      ymodMap = imDataMap(size(stackedImg),{gridy,gridx},gridYMod,'mask',cellMask);
   end

   bdfMap(find(bdfMap<0))   = 0;
   mcfMap(find(mcfMap<0))   = 0;
   adfMap(find(adfMap<0))   = 0;
   adcMap(find(adfMap<0))   = 0;
   spdMap(find(bdfMap<0))   = 0;
   mixfMap(find(mixfMap<0)) = 0;

   %We only show forces that are above a threshold.
   %We only show forces that are above a threshold.
   numInd = find(~isnan(bdfMap));
   if isempty(numInd)
      maxBDF = 0;
      avgBDF = 0;
   else
      maxBDF = max(bdfMap(numInd));
      avgBDF = mean(bdfMap(numInd));
   end
   smBDFInd = numInd(find(bdfMap(numInd)<avgBDF*smForceThreshold));
   bdfMap(smBDFInd) = NaN;

   numInd = find(~isnan(mcfMap));
   if isempty(numInd)
      maxMCF = 0;
      avgMCF = 0;
   else
      maxMCF = max(mcfMap(numInd));
      avgMCF = mean(mcfMap(numInd));
   end
   smMCFInd = numInd(find(mcfMap(numInd)<avgMCF*smForceThreshold));
   mcfMap(smMCFInd) = NaN;

   numInd = find(~isnan(adfMap));
   if isempty(numInd)
      maxADF = 0;
      avgADF = 0;
   else
      maxADF = max(adfMap(numInd));
      avgADF = mean(adfMap(numInd));
   end
   smADFInd = numInd(find(adfMap(numInd)<avgADF*smForceThreshold));
   adfMap(smADFInd) = NaN;
   adcMap(smADFInd) = NaN;

   %Identify the mixed zone.
   noMixZoneInd = find(isnan(mcfMap) | isnan(adfMap));
   mixfMap(noMixZoneInd) = NaN;
   
   %The displacements computed with the identified force and the identified
   % force at the data points. To be compared with 
   % 'dataU1' and 'dataU2'.
   %[dataUC1{jj} dataUC2{jj}] = postinterp(fem,'u1','u2', ...
   %   [dataPx{jj} dataPy{jj}].');
   %dataUC1{jj} = dataUC1{jj}.';
   %dataUC2{jj} = dataUC2{jj}.';

   %[dataBFx{jj} dataBFy{jj}] = postinterp(fem,'f1','f2', ...
   %   [dataPx{jj} dataPy{jj}].');
   %dataBFx{jj} = dataBFx{jj}.';
   %dataBFy{jj} = dataBFy{jj}.';

   %Calculate the residue of the forword computed displacements.
   %residueDisp{jj} = (dataUC1{jj}-dataU1{jj}).^2 + (dataUC2{jj}-dataU2{jj}).^2;

   %Calculate the residue of the regularized force.
   %residueBF{jj} = bfSigma*coef(:,jj).*coef(:,jj);
   %residueBF{jj} = bfSigma*(coefBFx(:,jj).*coefBFx(:,jj) + ...
   %   coefBFy(:,jj).*coefBFy(:,jj));

   %The displacement on the edge to be compared with 'edgeU1' etc.
   %for k = 1:numEdges
   %   [edgeUC1{jj,k} edgeUC2{jj,k}] = postinterp(fem,'u1','u2',edgeP{k});
   %end

   %Save the reconstructed forces.
   save(forceFieldFile,'forceField');
   save(iDispFieldFile,'iDispField');

   %Save the data map of reconstructed forces.
   indexStr = sprintf(imgIndexForm,imgIndexOfDTimePts(jj));
   bdfMapFile  = [reslDir filesep 'bdfMap' filesep 'bdfMap' indexStr '.mat'];
   mcfMapFile  = [reslDir filesep 'mcfMap' filesep 'mcfMap' indexStr '.mat'];
   adfMapFile  = [reslDir filesep 'adfMap' filesep 'adfMap' indexStr '.mat'];
   adcMapFile  = [reslDir filesep 'adcMap' filesep 'adcMap' indexStr '.mat'];
   mixfMapFile = [reslDir filesep 'mixfMap' filesep 'mixfMap' indexStr '.mat'];
   spdMapFile  = [reslDir filesep 'spdMap' filesep 'spdMap' indexStr '.mat'];
   ymodMapFile = [reslDir filesep 'ymodMap' filesep 'ymodMap' indexStr '.mat'];
   
   if ~exist('bdfMap','dir')
      [success msg msgID] = mkdir(reslDir,'bdfMap');
      if ~success
         error('Trouble making directory ''bdfMap''.');
      end
   end
   if ~exist('mcfMap','dir')
      [success msg msgID] = mkdir(reslDir,'mcfMap');
      if ~success
         error('Trouble making directory ''mcfMap''.');
      end
   end
   if ~exist('adfMap','dir')
      [success msg msgID] = mkdir(reslDir,'adfMap');
      if ~success
         error('Trouble making directory ''adfMap''.');
      end
   end
   if ~exist('adcMap','dir')
      [success msg msgID] = mkdir(reslDir,'adcMap');
      if ~success
         error('Trouble making directory ''adcMap''.');
      end
   end
   if ~exist('mixfMap','dir')
      [success msg msgID] = mkdir(reslDir,'mixfMap');
      if ~success
         error('Trouble making directory ''mixfMap''.');
      end
   end
   if ~exist('spdMap','dir')
      [success msg msgID] = mkdir(reslDir,'spdMap');
      if ~success
         error('Trouble making directory ''spdMap''.');
      end
   end
   if ~exist('ymodMap','dir')
      [success msg msgID] = mkdir(reslDir,'ymodMap');
      if ~success
         error('Trouble making directory ''ymodMap''.');
      end
   end
   save(bdfMapFile,'bdfMap');
   save(mcfMapFile,'mcfMap');
   save(adfMapFile,'adfMap');
   save(adcMapFile,'adcMap');
   save(mixfMapFile,'mixfMap');
   save(spdMapFile,'spdMap');
   save(ymodMapFile,'ymodMap');

   fprintf(1,'Done in %5.3f sec.\n',cputime-localStartTime);
end

%Save the identified body force calculated on the demonstration points.
%save([reslDir filesep 'bfId'],'bfDisplayPx','bfDisplayPy', ...
%   'recBF','recADFx','recADFy','recMCFx','recMCFy','recDispU','recStrain', ...
%   'recStress','residueDisp','residueBF', 'gridy','gridx','gridXin','gridYin', ...
%   'gridBFx','gridBFy', 'adfIndShow','mcfIndShow','mixIndShow', ...
%   'adfInd','mcfInd','mixInd','gridBDF','gridMCF','gridADF', ...
%   'gridADC','gridSpd','gridMixF','gridYMod');

%Save the computed displacement from the identified body force.
%save([reslDir filesep 'dispId'],'dataPx','dataPy','dataU1','dataU2', ...
%   'dataUC1','dataUC2', 'dataBFx','dataBFy');

fprintf(1,'Total time spent: %5.3f sec.\n',cputime-startTime);
