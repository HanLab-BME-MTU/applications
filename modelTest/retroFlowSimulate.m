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

if exist([simuDir filesep 'sDispField.mat'],'file') == 2
   answer = input(['A previously simulated flow field is found.\n' ...
      'Do you want to overwrite it? (y/n):'],'s');
   if strcmp(answer,'n')
      return;
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

   %fn.BndDispx{k} = 'bndDisp';
   %fn.BndDispy{k} = 'bndDisp';
   %fp.BndDispx{k} = {{'s'} {edgePPx{testTimeStep,k}}};
   %fp.BndDispy{k} = {{'s'} {edgePPy{testTimeStep,k}}};
end

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
   %Draw contraction region.
   answer = input('Do you want to draw a region of contraction? (y/n):','s');
   if strcmp(answer,'y')
      [BW,mcfROIXi,mcfROIYi] = roipoly;
   else
      mcfROIXi = [];
      mcfROIYi = [];
   end

   %Separate myosin contraction.
   mcfPx = forceField.p(:,1);
   mcfPy = forceField.p(:,2);
   mcfFx = forceField.mcf(:,1);
   mcfFy = forceField.mcf(:,2);
   nanInd = find(isnan(mcfFx) | isnan(mcfFy));
   mcfFx(nanInd) = 0;
   mcfFy(nanInd) = 0;

   %Calculate the DOF of the myosin dragging force in the fem stucture, 'fs.fem'.
   %coefMyo  = calFemDOFCoef(fs,[bfDisplayPx bfDisplayPy].',[mcfFx;mcfFy],1e-6);
   coefMyo  = calFemDOFCoef(fs.fem,[mcfPx mcfPy].',[mcfFx;mcfFy],1e-6);
   coefMyox = coefMyo(1,:).';
   coefMyoy = coefMyo(2,:).';

   fn.MyoDragFx = 'femBodyF';
   fn.MyoDragFy = 'femBodyF';
   fp.MyoDragFx = { {'x' 'y'} {fs.fem coefMyox mcfROIXi mcfROIYi} };
   fp.MyoDragFy = { {'x' 'y'} {fs.fem coefMyoy mcfROIXi mcfROIYi} };

   fn.BodyFx = 0;
   fn.BodyFy = 0;

   %Modify the body force.
   %fp.BodyFx = {{'x' 'y'} {fs coefBFx(:,testTimeStep) mcfROIXi mcfROIYi}};
   %fp.BodyFy = {{'x' 'y'} {fs coefBFy(:,testTimeStep) mcfROIXi mcfROIYi}};
end

bndfEdgNo = [];
if strcmp(forceToSimu,'bndOnly') == 1 || ...
   strcmp(forceToSimu,'adfAndBnd') == 1 || ...
   strcmp(forceToSimu,'mcfAndBnd') == 1 || ...
   strcmp(forceToSimu,'all') == 1
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
   options = elOptionsSet(options,'BCType',BCTypes);
end

if strcmp(forceToSimu,'adfOnly') == 1 || ...
   strcmp(forceToSimu,'adfAndBnd') == 1 || ...
   strcmp(forceToSimu,'mcfAndAdf') == 1 || ...
   strcmp(forceToSimu,'all') == 1

   %Load adhesion coefficients map.
   indexStr = sprintf(indexForm,imgIndexOfDTimePts(curDTimePt));
   adcMapFile = [reslDir filesep 'adcMap' filesep 'adcMap' indexStr '.mat'];
   load(adcMapFile);

   %Draw contraction region.
   answer = input('Do you want to draw a region of adhesion? (y/n):','s');
   if strcmp(answer,'y')
      [BW,adcROIXi,adcROIYi] = roipoly;
   else
      adcROIXi = [];
      adcROIYi = [];
   end

   imgHeight = size(adcMap,1);
   imgWidth  = size(adcMap,2);
   [pixX pixY] = meshgrid([1:imgWidth],[1:imgHeight]);
   in = inpolygon(pixX(:),pixY(:),adcROIXi,adcROIYi);
   iIn = find(in==1);
   adcMap(iIn) = 0;
   
   nanInd = find(isnan(adcMap));
   adcMap(nanInd) = 0;
   maxADC = max(adcMap(:));
   fn.VDragCoef = 'elImgParFun';
   fp.VDragCoef = { {'x' 'y'} {adcMap,0,maxADC,0,'Gaussian',2} };
end

fem = elModelUpdate(fem,'fn',fn,'fp',fp,'options',options);
fem = elasticSolve(fem,[]);

%Load the data points used to identify the force. We shall compute the 
% simulated displacements there.
[simU1 simU2]   = postinterp(fem,'u1','u2', forceField.p.');
[simBFx simBFy] = postinterp(fem,'f1','f2',[simDataPx simDataPy].');

sDispField.p   = forceField.p;
sDispField.v   = [simuU1;simuU2].';
sDispField.f   = [simuBFx;simuBFy].';
sDispField.fem = fem;

sDispField.bndfEdgNo = bndfEdgNo;
if ~isempty(bndfEdgNo)
   sDispField.bndF = forceField.bndF(bndfEdgNo);
end

save([simuDir filesep 'sDispField.mat'],'sDispField');


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
