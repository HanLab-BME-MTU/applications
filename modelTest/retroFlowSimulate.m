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

if isempty(simuDir) || ~isdir(simuDir)
   fprintf(1,'Please setup the simulation directory first by ''setupForceProj''.');
   return;
end

if exist([simuDir filesep 'simField.mat'],'file') == 2
   ans = input(['A previously simulated flow field is found.\n' ...
      'Do you want to load it? (y/n):'],'s');
   if strcmp(ans,'y')
      load([simuDir filesep 'simuField.mat']);
      return;
   end
end

%Load the 'fem' structure used to identify the force. This provides the same
% elasticity and boundary condtions.
load([mechDir filesep 'femId']);
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
load([reslDir filesep 'coefBFId']);
load([reslDir filesep 'bfId']);
load([reslDir filesep 'tfId']);

%Load the boundary displacements.
load([reslDir filesep 'edgeDisp']);
load([reslDir filesep 'dispId']);

fn.MyoDragFx = 0;
fn.MyoDragFy = 0;
fn.VDragCoef = 0;

fn.BodyFx = 0;
fn.BodyFy = 0;

%fn.BodyFx = 'femBodyF';
%fn.BodyFy = 'femBodyF';
%fp.BodyFx = {{'x' 'y'} {fs coefBFx(:,testTimeStep)}};
%fp.BodyFy = {{'x' 'y'} {fs coefBFy(:,testTimeStep)}};

if strcmp(forceToSimu,'mcfOnly') == 1 || ...
   strcmp(forceToSimu,'mcfAndBnd') == 1 || ...
   strcmp(forceToSimu,'mcfAndAdf') == 1 || ...
   strcmp(forceToSimu,'all') == 1
   %Draw contraction region.
   ans = input('Do you want to draw a region of contraction? (y/n):','s');
   if strcmp(ans,'y')
      [BW,mcfROIXi,mcfROIYi] = roipoly;
   else
      mcfROIXi = [];
      mcfROIYi = [];
   end

   %Separate myosin contraction.
   mcfFx = recMCFx;
   mcfFy = recMCFy;
   nanInd = find(isnan(recMCFx) | isnan(recMCFy));
   mcfFx(nanInd) = 0;
   mcfFy(nanInd) = 0;

   %Calculate the DOF of the myosin dragging force in the fem stucture, 'fs'.
   coefMyo  = calFemDOFCoef(fs,[bfDisplayPx bfDisplayPy].',[mcfFx;mcfFy],1e-6);
   coefMyox = coefMyo(1,:).';
   coefMyoy = coefMyo(2,:).';

   fn.MyoDragFx = 'femBodyF';
   fn.MyoDragFy = 'femBodyF';
   fp.MyoDragFx = { {'x' 'y'} {fs coefMyox mcfROIXi mcfROIYi} };
   fp.MyoDragFy = { {'x' 'y'} {fs coefMyoy mcfROIXi mcfROIYi} };

   fn.BodyFx = 0;
   fn.BodyFy = 0;

   %Modify the body force.
   %fp.BodyFx = {{'x' 'y'} {fs coefBFx(:,testTimeStep) mcfROIXi mcfROIYi}};
   %fp.BodyFy = {{'x' 'y'} {fs coefBFy(:,testTimeStep) mcfROIXi mcfROIYi}};
end

if strcmp(forceToSimu,'bndOnly') == 1 || ...
   strcmp(forceToSimu,'adfAndBnd') == 1 || ...
   strcmp(forceToSimu,'mcfAndBnd') == 1 || ...
   strcmp(forceToSimu,'all') == 1
   bndfEdgNo = input('Please choose the edges that you want to apply force to:\n');
   fprintf(1,['Specify forces:\n' ...
      '   1. Free edge.\n' '   2. Identified boundary force.\n']);
   for ii = 1:length(bndfEdgNo)
      k = bndfEdgNo(ii);
      BCTypes{k} = 'Neumann';

      ans = input(['Edge No. ' sprintf('%d',k) ': ']);
      if ans == 1
         fp.BndTracFx{k} = {{'s'} {0}};
         fp.BndTracFy{k} = {{'s'} {0}};
      elseif ans == 2
         fp.BndTracFx{k} = {{'s'} {spTFx{k}}};
         fp.BndTracFy{k} = {{'s'} {spTFy{k}}};
      end
   end
   options = elOptionsSet(options,'BCType',BCTypes);
end

if strcmp(forceToSimu,'adfOnly') == 1 || ...
   strcmp(forceToSimu,'adfAndBnd') == 1 || ...
   strcmp(forceToSimu,'mcfAndAdf') == 1 || ...
   strcmp(forceToSimu,'all') == 1

   %Load adhesion coefficients map.
   indexStr = sprintf(indexForm,imgIndexOfDTimePts(jj));
   adcMapFile = [reslDir filesep 'adcMap' filesep 'adcMap' indexStr '.mat'];
   load(adcMapFile);

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
simDataPx = dataPx{testTimeStep};
simDataPy = dataPy{testTimeStep};
[simU1 simU2] = postinterp(fem,'u1','u2', ...
   [simDataPx simDataPy].');
[simBFx simBFy] = postinterp(fem,'f1','f2',[simDataPx simDataPy].');

[simDispU1 simDispU2]     = postinterp(fem,'u1','u2',[bfDisplayPx bfDisplayPy].');
[simDispBFx simDispBFy]   = postinterp(fem,'f1','f2',[bfDisplayPx bfDisplayPy].');

for k = 1:numEdges
   %[simBndFx{k} simBndFy{k}] = postinterp(fem,'g1','g2',edgeS{k},'edim',1,'dom',k);
   simBndFx{k} = fnval(spTFx{k},edgeS{k});
   simBndFy{k} = fnval(spTFy{k},edgeS{k});
end

save([simuDir filesep 'simField'], 'fem','bfDisplayPx','bfDisplayPy', ...
   'simDispU1','simDispU2','simDispBFx','simDispBFy', ...
   'edgeP','simBndFx','simBndFy','simU1','simU2', ...
   'simDataPx','simDataPy','simBFx','simBFy');



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
