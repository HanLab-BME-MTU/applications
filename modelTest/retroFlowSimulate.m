%This script file simulates the retrograde flow with the force identified from
% real actin speckle movie. You have two options depending on the values of
% 'whatToTest':
%    'bodyForce' : The body force can be modified to see how it affects
%       the displacement field and the distorted displacement field 
%       can then be fed back into the force identification program to 
%       see if it can identify the modification we made on the force. 
%    'adhesion' : We can also use the adhesion intesity movies to add adhesion
%       site to see if they can be identified.

%Load the 'fem' structure used to identify the force. This provides the same
% elasticity and boundary condtions.
load([modelPath 'femId']);

%Load the boundary displacements.
load([resultPath 'edgeDisp']);
for k = 1:numEdges
   fn.BndDispx{k} = 'bndDisp';
   fn.BndDispy{k} = 'bndDisp';
   fp.BndDispx{k} = {{'s'} {edgePPx{testFrame,k}}};
   fp.BndDispy{k} = {{'s'} {edgePPy{testFrame,k}}};
end

if strcmp(whatToTest,'bodyForce') == 1
   %Load the identified body force.
   load([resultPath 'coefBFId']);
   fn.BodyFx = 'femBodyF';
   fn.BodyFy = 'femBodyF';
   fp.BodyFx = {{'x' 'y'} {fs coefBFx(:,testFrame)}};
   fp.BodyFy = {{'x' 'y'} {fs coefBFy(:,testFrame)}};

   %Modify the body force.
   fp.BodyFx = {{'x' 'y'} {fs coefBFx(:,testFrame)}};
   fp.BodyFy = {{'x' 'y'} {fs coefBFy(:,testFrame)}};
elseif strcmp(whatToTest,'adhesion') == 1
   load([resultPath 'myoDFId']);
   fn.BodyFx = 0;
   fn.BodyFy = 0;
   
   %Calculate the DOF of the myosin dragging force in the fem stucture, 'fs'.
   % Constructing the matrix whose columns are the values of the basis
   % function at the data points.
   dimFS  = flngdof(fs); 
   coefFS = zeros(dimFS,1);

   myoA = zeros(length(bfDisplayPx),dimFS);
   name = fs.dim;
   myoPx = bfDisplayPx;
   myoPy = bfDisplayPy;
   [is,pe] = postinterp(fs,[bfDisplayPx bfDisplayPy].');
   myoPx(pe) = [];
   myoPy(pe) = [];
   for k = 1:dimFS
      coefFS(k) = 1;
      myoA(:,k) = postinterp(fs,name,[myoPx myoPy].','u',coefFS);
      coefFS(k) = 0;
   end

   coefMyox = (myoA.'*myoA+1e-4*eye(dimFS))\(myoA.'*myoDFx);
   coefMyoy = (myoA.'*myoA+1e-4*eye(dimFS))\(myoA.'*myoDFy);

   fn.MyoDragFx = 'femBodyF';
   fn.MyoDragFy = 'femBodyF';
   fp.MyoDragFx = { {'x' 'y'} {fs coefMyox} };
   fp.MyoDragFy = { {'x' 'y'} {fs coefMyoy} };

   fn.VDragCoef = 'adhesionVC';
end

fem = elModelUpdate(fem,'fn',fn,'fp',fp);
fem = elasticSolve(fem,[]);

%Load the raw data points where we shall computed the simulated displacements.
load([resultPath 'dispField']);
simDataPx = rawDataPx{testFrame};
simDataPy = rawDataPy{testFrame};
[is,pe] = postinterp(fem,[simDataPx simDataPy].');
simDataPx(pe) = [];
simDataPy(pe) = [];
[simulU1 simulU2] = postinterp(fem,'u1','u2',[simDataPx simDataPy].');
%[simulBFx simulBFy] = postinterp(fem,'f1','f2',[dataPx dataPy].');
[simulBFx simulBFy] = postinterp(fem,'f1','f2',[bfDisplayPx bfDisplayPy].');

%Add noise to data.
%Noise level for testing simulation result.
noiseA = 0.0;
simulU1 = simulU1.*(1+noiseA*randn(size(simulU1)));
simulU2 = simulU2.*(1+noiseA*randn(size(simulU2)));

save([resultPath 'simField'],'simDataPx','simDataPy', ...
   'simulU1','simulU2','simulBFx','simulBFy');
