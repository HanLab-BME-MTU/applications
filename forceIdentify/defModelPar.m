%Set parameters for the finite element model.
meshType = 'rectMesh'; %Values: 'rectMesh' or 'femMesh'.
bFunType = 'femShape'; %Values: 'femShape' or 'BSpline'.
if strcmp(meshType,'rectMesh') == 1
   %Number of meshing points in the horizontal and vertical direction.
   recHin   = 60;
   recVin   = 40;
end

if strcmp(bFunType,'femShape')
   bfDomHin = 30;
   bfDomVin = 20;
end

%Specify elasticity
%fn.lambda  = 3;
%fn.mu      = 1;
%Add some variation to 'fn.YModul'.
%modulA = 2; %Average of modul.
%modulV = 0; %Magnitude of variation.
%varT   = 200; %Period of variation in pixels.
%fn.YModul = [modelPath 'YModulSine'];
%fp.YModul = { {'x'} {modulA,modulV,varT} };
fn.YModul = 2;
fn.PRatio = 0.3;
fp.YModul = [];
fp.PRatio = [];

%Specify the viscosity coefficient.
%fn.ViscCoef = {'viscCoef'};
fn.VDragCoef = 0.0; %0.0005;
fn.TimeStep  = 1;

modParam.meshType = meshType;
modParam.bFunType = bFunType;
modParam.recHin   = recHin;
modParam.recVin   = recVin;
modParam.bfDomHin = bfDomHin;
modParam.bfDomVin = bfDomVin;
modParam.fn       = fn;
modParam.fp       = fp;
