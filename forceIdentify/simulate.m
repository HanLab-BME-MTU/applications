clear all;

fprintf(1,'Start simulation ...\n');
startTime = cputime;

%load the geometry of the lamellipodium;
load data/simulGeom; %You get: 'geom', 'smplPGx' and 'smplPGy'.

%Get the starting and ending vertice number of each edge and the coordinate of
% each vertice.
[edgV vertP] = geominfo(geom,'out',{'se' 'mp'});

%Set up the edge label numbers for different parts of the boundary: the
% Leading Edge and the rest.
LEdgL = [8 7 6 4 2 1 3 5 14 20 21 23 27 31 38 42 44 45 43 41]; %Leading Edge.
FBndL = [9 10 15 18 25 30 32 36 40 35 33 28 24 22 16 12 11 13 ...
   17 19 26 29 34 37 35]; %The rest of the boundary, Free boundary.

%Group subdomains and boundaries.
ind = [1 2 3];
bndInd = ones(1,size(edgV,2));
bndInd(LEdgL) = 2;

options = elOptionsSet('EPType','YModulPRatio','BCType', ...
   {'Neumann' 'Neumann'});

%Specify elasticity
%fn.lambda  = 3;
%fn.mu      = 1;
fn.YModul = 2;
fn.PRatio = 0.3;

%Specify the viscosity coefficient.
viscA = 0.0005;
fn.VDragCoef = viscA;
%fn.VDragCoef = 'viscCoef';
%fp.VDragCoef = { {'x' 'y'} {viscA} };
fn.TimeStep  = 1;

%Specify the boundary traction force.
bndA = 5e-2;
fn.BndTracFx = {0 'bndTracFx'};
fn.BndTracFy = {0 'bndTracFy'};
fp.BndTracFx = { [] {{'x','nx'} {bndA}} };
fp.BndTracFy = { [] {{'x','ny'} {bndA}} };

%Specify boundary displacement.
fn.BndDispx = {[] []};
fn.BndDispy = {[] []};

%Myosin draging force.
myoA = 1e-2;
fn.MyoDragFx = {0 @myoDragFx 0};
fn.MyoDragFy = {0 'myoDragFy' 0};
fp.MyoDragFx = { [] {{'x','y'} {myoA}} [] };
fp.MyoDragFy = { [] {{'x','y'} {myoA}} [] };
%fn.MyoDragFx = {0 @myoDFx 0};
%fn.MyoDragFy = {0 'myoDFy' 0};
%fp.MyoDragFx = { [] {{'x','y'} {myoA}} [] };
%fp.MyoDragFy = { [] {{'x','y'} {0}} [] };

femS = elModelAssemble(geom,{'hmax',10},options,fn,fp,ind,bndInd);
femS = elasticSolve(femS,[]);

save data/simulData femS;
fprintf(1,'Done: %f sec.',cputime-startTime);
