%This script file simulates the retrograde flow with the force identified from
% real actin speckle movie. You have two options depending on the values of
% 'whatToTest':
%    'bodyForce' : The body force can be modified to see how it affects
%       the displacement field and the distorted displacement field 
%       can then be fed back into the force identification program to 
%       see if it can identify the modification we made on the force. 
%    'adhesion' : We can also use the adhesion intesity movies to add
%       viscosity dragging coefficient in the forward model to see if the
%       adhesion site can be identified.

%Load the 'fem' structure used to identify the force. This provides the same
% elasticity and boundary condtions.
load([modelPath 'femId']);

%Load the identified body force.
load([bfDataPath 'coefBFId']);

%Load the boundary displacements.
load([bfDataPath 'edgeDisp']);
for k = 1:numEdges
   fn.BndDispx{k} = 'bndDisp';
   fn.BndDispy{k} = 'bndDisp';
   fp.BndDispx{k} = {{'s'} {edgePPx{testTimeStep,k}}};
   fp.BndDispy{k} = {{'s'} {edgePPy{testTimeStep,k}}};
end

fn.MyoDragFx = 0;
fn.MyoDragFy = 0;
fn.VDragCoef = 0;

fn.BodyFx = 'femBodyF';
fn.BodyFy = 'femBodyF';
fp.BodyFx = {{'x' 'y'} {fs coefBFx(:,testTimeStep)}};
fp.BodyFy = {{'x' 'y'} {fs coefBFy(:,testTimeStep)}};

if strcmp(whatToTest,'bodyForce') == 1
   %Modify the body force.
   fp.BodyFx = {{'x' 'y'} {fs coefBFx(:,testTimeStep)}};
   fp.BodyFy = {{'x' 'y'} {fs coefBFy(:,testTimeStep)}};
elseif strcmp(whatToTest,'adhesion') == 1
   %Get the image of the adhesion intensity map that corresponds to the
   % 'testTimeStep'.
   img = imread(imgFile{startFrame+framesPerTStep*(testTimeStep-1)});

   %Set a threshold in terms of percentage where any adhesion intensity below
   % the threshold is set to be zero.
   adhThreshold = 0.1;

   %Compute the identified body force.
   fem = elModelUpdate(fem,'fn',fn,'fp',fp);
   fem = elasticSolve(fem,[]);

   load([bfDataPath 'bfId']);
   myoPx = [fs.mesh.p(1,:) bfDisplayPx.'];
   myoPy = [fs.mesh.p(2,:) bfDisplayPy.'];
   [is,pe] = postinterp(fs,[myoPx;myoPy]);
   myoPx(pe)  = [];
   myoPy(pe)  = [];

   [is,pe] = postinterp(fem,[myoPx;myoPy]);
   myoPx(pe)  = [];
   myoPy(pe)  = [];

   [recBFx, recBFy] = postinterp(fem,'f1','f2',[myoPx;myoPy]);
   [recDispU1, recDispU2] = postinterp(fem,'u1','u2',[myoPx;myoPy]);

   %Separate the body force into the viscosity dragging force and the
   % contraction force.
   recBFLen   = sqrt(recBFx.^2+recBFy.^2);
   recDispLen = sqrt(recDispU1.^2+recDispU2.^2);
   unitRecU1  = recDispU1./recDispLen;
   unitRecU2  = recDispU2./recDispLen;
   dotProdBFU = recBFx.*unitRecU1 + recBFy.*unitRecU2;

   % mcfInd : Indices where we consider the recovered body force to be mainly 
   %          myosin contraction.
   % adfInd : Indices where we consider the recovered body force to be mainly 
   %          adhesion.
   mcfInd = find((dotProdBFU-recBFLen*cos(mcfAngle))>=0);
   adfInd = find((-dotProdBFU-recBFLen*cos(adfAngle))>=0);
   myoDFx = zeros(size(recBFx));
   myoDFy = zeros(size(recBFy));
   myoDFx(mcfInd) = recBFx(mcfInd);
   myoDFy(mcfInd) = recBFy(mcfInd);

   %Set the myosin dragging force to be zero where the adhesion intensity is
   % higher than the threshold, 'adhThreshold'.
   % adhI      : Intensity of the adhesion at recovered body force points.
   % bigAdhInd : Indices where the adhesion is higher than the threshold
   adhI      = imInterp(img,[myoPy;myoPx].');
   bigAdhInd = find(adhI>adhThreshold);
   myoDFx(bigAdhInd) = 0;
   myoDFy(bigAdhInd) = 0;

   %Calculate the DOF of the myosin dragging force in the fem stucture, 'fs'.
   % Constructing the matrix whose columns are the values of the basis
   % function at the data points.
   % myoA : the matrix of the linear system.
   dimFS  = flngdof(fs); 
   coefFS = zeros(dimFS,1);

   name = fs.dim;

   myoA = zeros(length(myoPx),dimFS);
   for k = 1:dimFS
      coefFS(k) = 1;
      bspF = postinterp(fs,name,[myoPx;myoPy],'u',coefFS);
      myoA(:,k) = bspF.';
      coefFS(k) = 0;
   end

   coefMyox = (myoA.'*myoA+0.01*eye(dimFS))\(myoA.'*myoDFx.');
   coefMyoy = (myoA.'*myoA+0.01*eye(dimFS))\(myoA.'*myoDFy.');
   %coefMyox = myoA\(myoDFx.');
   %coefMyoy = myoA\(myoDFy.');

   fn.MyoDragFx = 'femBodyF';
   fn.MyoDragFy = 'femBodyF';
   fp.MyoDragFx = { {'x' 'y'} {fs coefMyox} };
   fp.MyoDragFy = { {'x' 'y'} {fs coefMyoy} };

   fn.BodyFx = 0;
   fn.BodyFy = 0;
   
   %Compute the maximum viscosity dragging coefficient.
   smDispLen = max(recDispLen)*smDispThreshold;
   smDispInd = find(recDispLen<smDispLen);
   recDispLen(smDispInd) = smDispLen;
   maxViscDC = max(recBFLen(adfInd)./recDispLen(adfInd));

   fn.VDragCoef = 'elImgParFun';
   fp.VDragCoef = { {'x' 'y'} {img,0,maxViscDC,adhThreshold} };
end

fem = elModelUpdate(fem,'fn',fn,'fp',fp);
fem = elasticSolve(fem,[]);

%Load the data points used to identify the force. We shall compute the 
% simulated displacements there.
load([bfDataPath 'dispId']);
simDataPx = dataPx{testTimeStep};
simDataPy = dataPy{testTimeStep};
[simulU1 simulU2] = postinterp(fem,'u1','u2', ...
   [simDataPx simDataPy].');
%[simulBFx simulBFy] = postinterp(fem,'f1','f2',[dataPx dataPy].');
[simMyoDFx simMyoDFy] = postinterp(fem,'f1','f2',[bfDisplayPx bfDisplayPy].');
[simDispU1 simDispU2] = postinterp(fem,'u1','u2',[bfDisplayPx bfDisplayPy].');

save([resultPath 'simField'], 'simDataPx','simDataPy','simulU1','simulU2', ...
   'simDispU1','simDispU2','simMyoDFx','simMyoDFy');
