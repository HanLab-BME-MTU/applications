%This script file simulates the retrograde flow with the force identified from
% real actin speckle movie. You have two options depending on the values of
% 'forceToSimu':
%    'bodyForce' : The body force can be modified to see how it affects
%       the displacement field and the distorted displacement field 
%       can then be fed back into the force identification program to 
%       see if it can identify the modification we made on the force. 
%    'adhesion' : We can also use the adhesion intesity movies to add
%       viscosity dragging coefficient in the forward model to see if the
%       adhesion site can be identified.
%
% AUTHOR: Lin Ji
% DATE  : July 17, 2006

if isempty(simuDir) || ~isdir(simuDir)
   fprintf(1,'Please setup the simulation directory first by ''setupForceProj''.');
   return;
end

if ~isdir([simuDir filesep 'rawDispField'])
   success = mkdir(simuDir,'rawDispField');
   if ~success
      error('Trouble making directory.');
   end
end

if ~isdir([simuDir filesep 'forceField'])
   success = mkdir(simuDir,'forceField');
   if ~success
      error('Trouble making directory.');
   end
end

numNoiseLevels = length(simRelNoiseLevel);
simIndexForm  = sprintf('%%.%dd',max(length(num2str(numNoiseLevels)),2));

sDispFieldDir = [simuDir filesep 'rawDispField'];
sForceFieldDir = [simuDir filesep 'forceField'];

if exist([sDispFieldDir filesep 'rawDispField.mat'],'file') == 2
   answer = input(['A previously simulated flow field is found.\n' ...
      'Do you want to redo the simulation? (y/n):'],'s');
   if strcmp(answer,'n')
      s = load([sDispFieldDir filesep 'rawDispField.mat']); 
      rawDispField = s.rawDispField;
      addSimNoise(rawDispField,sDispFieldDir,simAbsNoiseLevel,simRelNoiseLevel);
      return;
   end
end

if exist('rawDispField','var')
   clear rawDispField;
end

if exist('ROI','var')
   clear ROI;
end

ROI.bndfEdgNo = [];
ROI.mcfXi  = [];
ROI.mcfYi  = [];
ROI.adcXi  = [];
ROI.adcYi  = [];

if exist([simuDir filesep 'ROI.mat'],'file') == 2
   answer = input(['A previously drawn ROI is found.\n' ...
      'Do you want to load it? (y/n):'],'s');
   if strcmp(answer,'y')
      s = load([simuDir filesep 'ROI.mat']); 
      ROI = s.ROI;
   end
end

imgIndex = imgIndexOfDTimePts(curDTimePt);

%Load the 'fem' structure used to identify the force. This provides the same
% elasticity and boundary condtions.
if strcmp(isFieldBndFixed,'yes')
   femModelFile = [femModelDir filesep 'femModel' sprintf(imgIndexForm,0) '.mat'];
else
   [femModelFile,femModelImgIndex] = getFemModelFile(femModelDir,imgIndex,imgIndexForm);
end
s = load(femModelFile);
femModel = s.femModel;
fem      = femModel.fem;
options  = femModel.options;
fs       = femModel.fs;
numEdges = femModel.numEdges;

fn = fem.fn;
fp = fem.fp;

fn.BndDispx  = cell(1,numEdges);
fn.BndDispy  = cell(1,numEdges);
fp.BndDispx  = cell(1,numEdges);
fp.BndDispy  = cell(1,numEdges);
fn.BndTracFx = cell(1,numEdges);
fn.BndTracFy = cell(1,numEdges);
fp.BndTracFx = cell(1,numEdges);
fp.BndTracFy = cell(1,numEdges);

for k = 1:numEdges
   fn.BndDispx{k} = 0;
   fn.BndDispy{k} = 0;

   %Set up link to boundary traction force function.
   fn.BndTracFx{k} = 'spBndTF';
   fn.BndTracFy{k} = 'spBndTF';
   fp.BndTracFx{k} = {{'s'} {0}};
   fp.BndTracFy{k} = {{'s'} {0}};

   %Specify the type of boundary conditions.
   BCTypes{k} = 'Dirichlet';
end

%Load the displacement field used to reconstruct the force we are going to modify for simulation.
% The flow field of the 'rawDispField' structure, v, will be replaced by the simulated flow. And,
% the smoothed flow field, sv will be removed. Some new fields will be added to the structure
% especially the simulation force field, f.
rawDispFieldFileName = ['rawDispField' sprintf(imgIndexForm,imgIndex) '.mat'];
rawDispFieldFile     = [rawDispFieldDir filesep rawDispFieldFileName];
s = load(rawDispFieldFile);
rawDispField = s.rawDispField;

%Load the identified body force and boundary force.
forceFieldFile = [forceFieldDir filesep 'forceField' ...
   sprintf(imgIndexForm,imgIndex) '.mat'];
s = load(forceFieldFile);
forceField = s.forceField;

fn.MyoDragFx = 0;
fn.MyoDragFy = 0;
fn.VDragCoef = 0;

fn.BodyFx = 0;
fn.BodyFy = 0;

if strcmp(forceToSimu,'mcfOnly') == 1 || ...
   strcmp(forceToSimu,'mcfAndBnd') == 1 || ...
   strcmp(forceToSimu,'mcfAndAdf') == 1 || ...
   strcmp(forceToSimu,'all') == 1

   cellImg = imread(imgFileList{imgChannel}{imgIndex});
   figure; imshow(cellImg,[]); hold on;
   quiver(forceField.p(:,1),forceField.p(:,2),forceField.f(:,1)*bfScale, ...
      forceField.f(:,2)*bfScale,0,'r');

   if ~isempty(ROI.mcfXi)
      for kk = 1:length(ROI.mcfXi)
         plot(ROI.mcfXi{kk},ROI.mcfYi{kk},'w-.');
      end
   end

   %Draw contraction region.
   answer = input('Do you want to redraw regions of contraction? (y/n):','s');
   if strcmp(answer,'y')
      ROI.mcfXi = [];
      ROI.mcfYi = [];
      numMCFROIs = input('How many regions do you want to draw?: ');
      if numMCFROIs> 0
         for kk = 1:numMCFROIs
            [BW,mcfXi,mcfYi] = roipoly;
            ROI.mcfXi = [ROI.mcfXi {mcfXi}];
            ROI.mcfYi = [ROI.mcfYi {mcfYi}];
            plot(mcfXi,mcfYi,'g-.');
         end
      end
      save([simuDir filesep 'ROI.mat'],'ROI');
   end

   %Separate myosin contraction.
   mcfPx = forceField.p(:,1);
   mcfPy = forceField.p(:,2);
   mcfFx = simMCFMagFactor*forceField.mcf(:,1);
   mcfFy = simMCFMagFactor*forceField.mcf(:,2);
   nanInd = find(isnan(mcfFx) | isnan(mcfFy));
   mcfPx(nanInd) = [];
   mcfPy(nanInd) = [];
   mcfFx(nanInd) = [];
   mcfFy(nanInd) = [];

   %Calculate the DOF of the myosin dragging force in the fem stucture, 'fs.fem'.
   %coefMyo  = calFemDOFCoef(fs,[bfDisplayPx bfDisplayPy].',[mcfFx;mcfFy],1e-6);
   coefMyo  = calFemDOFCoef(fs.fem,[mcfPx mcfPy].',[mcfFx mcfFy].',1e-6);
   coefMyox = coefMyo(1,:).';
   coefMyoy = coefMyo(2,:).';

   fn.MyoDragFx = 'femBodyF';
   fn.MyoDragFy = 'femBodyF';
   fp.MyoDragFx = { {'x' 'y'} {fs.fem coefMyox ROI.mcfXi ROI.mcfYi} };
   fp.MyoDragFy = { {'x' 'y'} {fs.fem coefMyoy ROI.mcfXi ROI.mcfYi} };

   fn.BodyFx = 0;
   fn.BodyFy = 0;

   %Modify the body force.
   %fp.BodyFx = {{'x' 'y'} {fs coefBFx(:,testTimeStep) mcfXi mcfYi}};
   %fp.BodyFy = {{'x' 'y'} {fs coefBFy(:,testTimeStep) mcfXi mcfYi}};
end

if strcmp(forceToSimu,'bndOnly') == 1 || ...
   strcmp(forceToSimu,'adfAndBnd') == 1 || ...
   strcmp(forceToSimu,'mcfAndBnd') == 1 || ...
   strcmp(forceToSimu,'all') == 1

   if ~isempty(ROI.bndfEdgNo)
      fprintf(1,'Previously chosen edges: %d.\n',ROI.bndfEdgNo);
   end

   bndfEdgNo = input('Please choose the edges that you want to apply force to:\n');
   fprintf(1,['Specify forces:\n' ...
      '   1. Free edge.\n' '   2. Identified boundary force.\n']);
   for ii = 1:length(bndfEdgNo)
      k = bndfEdgNo(ii);
      bndF = forceField.bndF(k);

      BCTypes{k} = 'Neumann';

      answer = input(['Edge No. ' sprintf('%d',k) ': ']);
      if answer == 1
         fp.BndTracFx{k} = {{'s'} {0}};
         fp.BndTracFy{k} = {{'s'} {0}};
      elseif answer == 2
         fp.BndTracFx{k} = {{'s'} {bndF.spx}};
         fp.BndTracFy{k} = {{'s'} {bndF.spy}};
      end
   end

   ROI.bndfEdgNo = bndfEdgNo;
   save([simuDir filesep 'ROI.mat'],'ROI');

   options = elOptionsSet(options,'BCType',BCTypes);
end

if strcmp(forceToSimu,'adfOnly') == 1 || ...
   strcmp(forceToSimu,'adfAndBnd') == 1 || ...
   strcmp(forceToSimu,'mcfAndAdf') == 1 || ...
   strcmp(forceToSimu,'all') == 1

   %Load adhesion coefficients map.
   indexStr = sprintf(imgIndexForm,imgIndex);
   adcMapFile = [reslDir filesep 'adcMap' filesep 'adcMap' indexStr '.mat'];
   load(adcMapFile);

   if ~isempty(ROI.adcXi)
      for kk = 1:length(ROI.adcXi)
         plot(ROI.adcXi{kk},ROI.adcYi{kk},'w-.');
      end
   end

   %Draw adhesion region.
   answer = input('Do you want to redraw regions of adhesion? (y/n):','s');
   if strcmp(answer,'y')
      ROI.adcXi = [];
      ROI.adcYi = [];
      numADCROIs = input('How many regions do you want to draw?: ');
      if numADCROIs> 0
         for kk = 1:numADCROIs
            [BW,adcXi,adcYi] = roipoly;
            ROI.adcXi = [ROI.adcXi {adcXi}];
            ROI.adcYi = [ROI.adcYi {adcYi}];
            plot(adcXi,adcYi,'g-.');
         end
      end
      save([simuDir filesep 'ROI.mat'],'ROI');
   end

   imgHeight = size(adcMap,1);
   imgWidth  = size(adcMap,2);
   [pixX pixY] = meshgrid([1:imgWidth],[1:imgHeight]);
   if ~isempty(ROI.adcXi)
      iOut = 1:length(pixX(:));
      for kk = 1:length(ROI.adcXi)
         in  = inpolygon(pixX(:),pixY(:),ROI.adcXi{kk},ROI.adcYi{kk});
         iIn = find(in==1);
         iOut(iIn) = 0;
      end
      iOut(find(iOut==0)) = [];
      adcMap(iOut) = 0;
   end
   
   nanInd = find(isnan(adcMap));
   adcMap(nanInd) = 0;
   adcMap = adcMap*simADCMagFactor;
   maxADC = max(adcMap(:));
   fn.VDragCoef = 'elImgParFun';
   fp.VDragCoef = { {'x' 'y'} {adcMap,0,maxADC,0,'Gaussian',2} };
end

fem = elModelUpdate(fem,'fn',fn,'fp',fp,'options',options);
fem = elasticSolve(fem,[]);

%Load the data points used to identify the force. We shall compute the 
% simulated displacements there.
%First, get data points that are inside the mesh.
[is,pe] = postinterp(fem,rawDispField.p.');
rawDispField.p(pe,:) = [];

[simU1 simU2]     = postinterp(fem,'u1','u2', rawDispField.p.');
[simMCFx simMCFy] = postinterp(fem,'f1','f2', rawDispField.p.');

%Calculate adhesion resistance force.
VDragCoef = postinterp(fem,'VDragCoef',rawDispField.p.');
simADFx   = -VDragCoef.*simU1/fn.TimeStep;
simADFy   = -VDragCoef.*simU2/fn.TimeStep;

simU   = [simU1;simU2].';
simSpd = sqrt(sum(simU.^2,2));
medSpd = median(simSpd);

rawDispField.v      = simU;
rawDispField.medSpd = medSpd;
rawDispField.mcf    = [simMCFx;simMCFy].';
rawDispField.adf    = [simADFx;simADFy].';
rawDispField.f      = rawDispField.mcf + rawDispField.adf;
rawDispField.fem    = fem;

%Calculate edge displacement.
numEdges = femModel.numEdges;
edge     = femModel.edge;

for kk = 1:numEdges
   %[simBndU1 simBndU2] = postinterp(fem,'u1','u2',bndS,'edim',1,'dom',kk);
   [simBndU1 simBndU2] = postinterp(fem,'u1','u2',forceField.bndF(kk).p,'Ext',1);
   numInd = find(~isnan(simBndU1) & ~isnan(simBndU2));
   edgD(kk).U1 = simBndU1(numInd).';
   edgD(kk).U2 = simBndU2(numInd).';
   edgD(kk).p  = forceField.bndF(kk).p(:,numInd);

   %Create spline interpolation of the edge displacement using arclength
   % parameter stored in 'msh.e(3:4,:)'.
   edgD(kk).s    = forceField.bndF(kk).s(numInd);
   edgD(kk).sKnt = augknt(edgD(kk).s,2);
   edgD(kk).ppU1 = spapi(edgD(kk).sKnt,edgD(kk).s,edgD(kk).U1.');
   edgD(kk).ppU2 = spapi(edgD(kk).sKnt,edgD(kk).s,edgD(kk).U2.');
end
rawDispField.edgD = edgD;

if ~isempty(ROI.bndfEdgNo)
   rawDispField.bndF = forceField.bndF(bndfEdgNo);
else
   rawDispField.bndF = [];
end

%Save the flow field without noise.
rawDispField.absNoiseLevel = 0;
rawDispField.relNoiseLevel = 0;
sDispFieldFile = [sDispFieldDir filesep 'rawDispField.mat'];
save(sDispFieldFile,'rawDispField');

addSimNoise(rawDispField,sDispFieldDir,simAbsNoiseLevel,simRelNoiseLevel);

%Calculate forces on grid points so that color map of the true force can be produced.
iDispFieldFileName = ['iDispField' sprintf(imgIndexForm,imgIndex) '.mat'];
iDispFieldFile     = [iDispFieldDir filesep iDispFieldFileName];
s = load(iDispFieldFile);
iDispField = s.iDispField;

gridx = iDispField.gridx;
gridy = iDispField.gridy;
gridX = iDispField.gridX;
gridY = iDispField.gridY;

[is,pe]    = postinterp(fem,[gridX gridY].');
gridIn     = 1:length(gridX);
gridIn(pe) = [];
gridXin    = gridX(gridIn);
gridYin    = gridY(gridIn);

gridMCF = NaN*ones(length(gridy)*length(gridx),1);
gridADF = NaN*ones(length(gridy)*length(gridx),1);
gridSpd = NaN*ones(length(gridy)*length(gridx),1);

[gridU1,gridU2]     = postinterp(fem,'u1','u2',[gridXin gridYin].');
[gridMCFx,gridMCFy] = postinterp(fem,'f1','f2',[gridXin gridYin].');

VDragCoef = postinterp(fem,'VDragCoef',[gridXin gridYin].');
gridADFx = -VDragCoef.*gridU1/fn.TimeStep;
gridADFy = -VDragCoef.*gridU2/fn.TimeStep;

gridMCF(gridIn) = sqrt(gridMCFx.^2 + gridMCFy.^2);
gridADF(gridIn) = sqrt(gridADFx.^2 + gridADFy.^2);
gridSpd(gridIn) = sqrt(gridU1.^2 + gridU2.^2);

gridMCF = reshape(gridMCF,length(gridy),length(gridx));
gridADF = reshape(gridADF,length(gridy),length(gridx));
gridSpd = reshape(gridSpd,length(gridy),length(gridx));
gridBDF = gridMCF + gridADF;

%Get the mask of cell from 'edgeDir'. It will be used in creating data map.
cellMask = [];
if isdir(edgeDir)
   maskDir = [edgeDir filesep 'cell_mask'];
   if isdir(maskDir)
      [maskFileList maskIndex] = getNamedFiles(maskDir,'mask_');

      thisMaskIndex = find(maskIndex==imgIndex);
      if ~isempty(thisMaskIndex)
         cellMaskFile = maskFileList{thisMaskIndex(1)};
         cellMask = imread([maskDir filesep cellMaskFile]);
      end
   end
end

if isempty(cellMask)
   bdfMap = imDataMap(size(cellImg),{gridy,gridx},gridBDF,'bnd',[bfDomPGx bfDomPGy]);
else
   bdfMap = imDataMap(size(cellImg),{gridy,gridx},gridBDF,'mask',cellMask);
end

if isempty(cellMask)
   mcfMap  = imDataMap(size(cellImg),{gridy,gridx},gridMCF,'bnd',[bfDomPGx bfDomPGy]);
else
   mcfMap  = imDataMap(size(cellImg),{gridy,gridx},gridMCF,'mask',cellMask);
end

if isempty(cellMask)
   adfMap  = imDataMap(size(cellImg),{gridy,gridx},gridADF,'bnd',[bfDomPGx bfDomPGy]);
else
   adfMap  = imDataMap(size(cellImg),{gridy,gridx},gridADF,'mask',cellMask);
end

if isempty(cellMask)
   spdMap  = imDataMap(size(cellImg),{gridy,gridx},gridSpd,'bnd',[bfDomPGx bfDomPGy]);
else
   spdMap  = imDataMap(size(cellImg),{gridy,gridx},gridSpd,'mask',cellMask);
end

bdfMap(find(bdfMap<0)) = 0;
mcfMap(find(mcfMap<0)) = 0;
adfMap(find(adfMap<0)) = 0;
spdMap(find(spdMap<0)) = 0;

mixfMap = NaN*ones(size(adfMap));

%Save the data map of true forces.
indexStr = sprintf(simIndexForm,[]);
bdfMapFile  = [simuDir filesep 'bdfMap' filesep 'bdfMap' indexStr '.mat'];
mcfMapFile  = [simuDir filesep 'mcfMap' filesep 'mcfMap' indexStr '.mat'];
adfMapFile  = [simuDir filesep 'adfMap' filesep 'adfMap' indexStr '.mat'];
spdMapFile  = [simuDir filesep 'spdMap' filesep 'spdMap' indexStr '.mat'];
mixfMapFile = [simuDir filesep 'mixfMap' filesep 'mixfMap' indexStr '.mat'];

if ~exist('bdfMap','dir')
   [success msg msgID] = mkdir(simuDir,'bdfMap');
   if ~success
      error('Trouble making directory ''bdfMap''.');
   end
end
if ~exist('mcfMap','dir')
   [success msg msgID] = mkdir(simuDir,'mcfMap');
   if ~success
      error('Trouble making directory ''mcfMap''.');
   end
end
if ~exist('adfMap','dir')
   [success msg msgID] = mkdir(simuDir,'adfMap');
   if ~success
      error('Trouble making directory ''adfMap''.');
   end
end
if ~exist('spdMap','dir')
   [success msg msgID] = mkdir(simuDir,'spdMap');
   if ~success
      error('Trouble making directory ''spdMap''.');
   end
end
if ~exist('mixfMap','dir')
   [success msg msgID] = mkdir(simuDir,'mixfMap');
   if ~success
      error('Trouble making directory ''mixfMap''.');
   end
end

save(bdfMapFile,'bdfMap');
save(mcfMapFile,'mcfMap');
save(adfMapFile,'adfMap');
save(spdMapFile,'spdMap');
save(mixfMapFile,'mixfMap');

sForceField.p       = forceField.p;
sForceField.f       = rawDispField.f;
sForceField.bndF    = rawDispField.bndF;
sForceField.gridIn  = forceField.gridIn;
sForceField.mcf     = rawDispField.mcf;
sForceField.adf     = rawDispField.adf;
sForceField.gridMCF = gridMCF;
sForceField.gridBDF = gridBDF;
sForceField.gridADF = gridADF;
sForceField.gridSpd = gridSpd;

forceField = sForceField;
sForceFieldFile = [sForceFieldDir filesep 'forceField.mat'];
save(sForceFieldFile,'forceField');

sDispField.p  = iDispField.p;
sDispField.rv = rawDispField.v;
sDispField.gridx = iDispField.gridx;
sDispField.gridy = iDispField.gridy;
sDispField.gridX = iDispField.gridX;
sDispField.gridY = iDispField.gridY;

iDispField = sDispField;
sDispFieldFile = [simuDir filesep 'iDispField' filesep 'iDispField.mat'];
save(sDispFieldFile,'iDispField');

%%%% Some old code %%%%%%%%
   %if strcmp(forceToSimu,'adfOnly') == 1 || ...
   %   strcmp(forceToSimu,'adfAndBnd') == 1 || ...
%   strcmp(forceToSimu,'mcfAndAdf') == 1 || ...
%   strcmp(forceToSimu,'all') == 1
%   %Get the image of the adhesion intensity map that corresponds to the
%   % 'testTimeStep'.
%   img = imread(imgFile{startFrame+framesPerTStep*(testTimeStep-1)});
%
%   %Set a threshold in terms of percentage where any adhesion intensity below
%   % the threshold is set to be zero.
%   adhThreshold = 0.1;
%
%   %Compute the identified body force.
%   fem = elModelUpdate(fem,'fn',fn,'fp',fp);
%   fem = elasticSolve(fem,[]);
%
%   load([bfDataPath 'bfId']);
%   myoPx = [fs.mesh.p(1,:) bfDisplayPx.'];
%   myoPy = [fs.mesh.p(2,:) bfDisplayPy.'];
%   [is,pe] = postinterp(fs,[myoPx;myoPy]);
%   myoPx(pe)  = [];
%   myoPy(pe)  = [];
%
%   [is,pe] = postinterp(fem,[myoPx;myoPy]);
%   myoPx(pe)  = [];
%   myoPy(pe)  = [];
%
%   [recBFx, recBFy] = postinterp(fem,'f1','f2',[myoPx;myoPy]);
%   [recDispU1, recDispU2] = postinterp(fem,'u1','u2',[myoPx;myoPy]);
%
%   %Separate the body force into the viscosity dragging force and the
%   % contraction force.
%   recBFLen   = sqrt(recBFx.^2+recBFy.^2);
%   recDispLen = sqrt(recDispU1.^2+recDispU2.^2);
%   unitRecU1  = recDispU1./recDispLen;
%   unitRecU2  = recDispU2./recDispLen;
%   dotProdBFU = recBFx.*unitRecU1 + recBFy.*unitRecU2;
%
%   % mcfInd : Indices where we consider the recovered body force to be mainly 
%   %          myosin contraction.
%   % adfInd : Indices where we consider the recovered body force to be mainly 
%   %          adhesion.
%   mcfInd = find((dotProdBFU-recBFLen*cos(mcfAngle))>=0);
%   adfInd = find((-dotProdBFU-recBFLen*cos(adfAngle))>=0);
%   myoDFx = zeros(size(recBFx));
%   myoDFy = zeros(size(recBFy));
%   myoDFx(mcfInd) = recBFx(mcfInd);
%   myoDFy(mcfInd) = recBFy(mcfInd);
%
%   if exist('bfConstraint') == 1 & strcmp(bfConstraint,'adhLocation') == 1
%      %Set the myosin dragging force to be zero at the adhesion sites
%      % segmented in 'adhLabel' or defined by 'adhGeom'.
%      if strcmp(adhGeomType,'node') == 1
%         load([modelPath 'adhSeg']);
%
%         %At adhesion pixels, the segmented and label image is non-zero.
%         adhInd = find(adhLabel(ceil(myoPy) + ...
%            size(adhLabel,1)*(ceil(myoPx)-1))~=0);
%         myoDFx(adhInd) = 0;
%         myoDFy(adhInd) = 0;
%
%         %Set the adhesion intensity to be zero outside the adhesion site.
%         img(find(adhLabel==0)) = 0;
%      elseif strcmp(adhGeomType,'convexHull') == 1
%         load([modelPath 'adhGeom']);
%
%         % If the adhesion site is defined as a geometrical object, we first 
%         % build a 'fem' structure for the adhesion geometry.
%         adhFem.geom  = adhGeom;
%         adhFem.mesh  = meshinit(adhFem);
%         adhFem.xmesh = meshextend(adhFem);
%
%         [is,pe] = postinterp(adhFem,[myoPx;myoPy]);
%         %'pe' : The indices of '[myoPx;myoPy]' that are outside 'adhGeom'.
%         %'adhInd' : The indices of '[myoPx;myoPy]' that are inside 'adhGeom'.
%         adhInd = [1:length(myoPx)];
%         adhInd(pe) = [];
%         myoDFx(adhInd) = 0;
%         myoDFy(adhInd) = 0;
%
%         %Set the adhesion intensity to be zero outside 'adhGeom'.
%         %Number of pixels in the horizontal direction
%         numPixelsH = size(img,2); 
%
%         %Number of pixels in the vertical direction
%         numPixelsV = size(img,1); 
%
%         [pixelH,pixelV] = meshgrid([1:numPixelsH],[1:numPixelsV]);
%
%         pixelH = reshape(pixelH,1,length(pixelH(:)));
%         pixelV = reshape(pixelV,1,length(pixelV(:)));
%
%         [is,pe] = postinterp(adhFem,[pixelH;pixelV]);
%         img(pe) = 0;
%      end
%   else
%      %Set the myosin dragging force to be zero where the adhesion intensity is
%      % higher than the threshold, 'adhThreshold'.
%      % adhI      : Intensity of the adhesion at recovered body force points.
%      % bigAdhInd : Indices where the adhesion is higher than the threshold
%      adhI      = imInterp(img,[myoPy;myoPx].');
%      bigAdhInd = find(adhI>adhThreshold);
%      myoDFx(bigAdhInd) = 0;
%      myoDFy(bigAdhInd) = 0;
%   end
%
%   %Calculate the DOF of the myosin dragging force in the fem stucture, 'fs'.
%   coefMyo  = calFemDOFCoef(fs,[myoPx;myoPy],[myoDFx;myoDFy],0.01);
%   coefMyox = coefMyo(1,:).';
%   coefMyoy = coefMyo(2,:).';
%
%   fn.MyoDragFx = 'femBodyF';
%   fn.MyoDragFy = 'femBodyF';
%   fp.MyoDragFx = { {'x' 'y'} {fs coefMyox} };
%   fp.MyoDragFy = { {'x' 'y'} {fs coefMyoy} };
%
%   fn.BodyFx = 0;
%   fn.BodyFy = 0;
%   
%   %Compute the maximum viscosity dragging coefficient.
%   smDispLen = max(recDispLen)*smDispThreshold;
%   smDispInd = find(recDispLen<smDispLen);
%   recDispLen(smDispInd) = smDispLen;
%   maxViscDC = max(recBFLen(adfInd)./recDispLen(adfInd));
%
%   fn.VDragCoef = 'elImgParFun';
%   fp.VDragCoef = { {'x' 'y'} {img,0,maxViscDC,adhThreshold} };
%end
