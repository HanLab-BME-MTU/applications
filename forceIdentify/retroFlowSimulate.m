%This script file simulates the retrograde flow with the force identified from
% real actin speckle movie. The force can be modified to see how it affects
% the displacement field and the distorted displacement field can then be fed
% back into the force identification program to see if it can identify the
% modification we made on the force.

%Load the 'fem' structure used to identify the force. This provides the same
% elasticity and boundary condtions.
load([modelPath 'femId']);

%Load the identified body force.
load([resultPath 'bfId']);
fn.BodyFx = 'spMyoDFx';
fn.BodyFy = 'spMyoDFy';
fp.BodyFx = {{'x' 'y'} {fs coefBF(1:end/2)}};
fp.BodyFy = {{'x' 'y'} {fs coefBF(end/2+1:end)}};

%Load the boundary displacements.
load([resultPath 'edgeDisp']);
for k = 1:numEdges
   fn.BndDispx{k} = 'bndDisp';
   fn.BndDispy{k} = 'bndDisp';
   fp.BndDispx{k} = {{'s'} {edgePPx{k}}};
   fp.BndDispy{k} = {{'s'} {edgePPy{k}}};
end

%Modify the body force.
fp.BodyFx = {{'x' 'y'} {fs coefBF(1:end/2)/2}};
fp.BodyFy = {{'x' 'y'} {fs coefBF(end/2+1:end)}};

fem = elModelUpdate(fem,'fn',fn,'fp',fp);
fem = elasticSolve(fem,[]);

%Load the raw data points where we shall computed the simulated displacements.
load([resultPath 'dispField']);
dataPx = rawDataPx;
dataPy = rawDataPy;
[is,pe] = postinterp(fem,[rawDataPx rawDataPy].');
dataPx(pe) = [];
dataPy(pe) = [];
[simulU1 simulU2] = postinterp(fem,'u1','u2',[dataPx dataPy].');
[simulBFx simulBFy] = postinterp(fem,'f1','f2',[dataPx dataPy].');

%Add noise to data.
%Noise level for testing simulation result.
noiseA = 0.0;
simulU1 = simulU1.*(1+noiseA*randn(size(simulU1)));
simulU2 = simulU2.*(1+noiseA*randn(size(simulU2)));

save([resultPath 'simField'],'dataPx','dataPy','simulU1','simulU2');
