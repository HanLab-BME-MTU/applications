clear all;

%Construct the geometry first.
%Load the data points of the boundary.
%load bndDataP;
%c1 = geomspline(p1,'splinemethod','centripetal','closed','off');
%c2 = geomspline(p2,'splinemethod','foley','closed','off');
%c3 = geomspline(p3,'splinemethod','chordlength','closed','off');
%c  = geomcoerce('solid',{c1,c2,c3});

c1 = curve2([-1 0 1],[0 2/3 0]);
c2 = curve2([-1 0 1],[0 1/6 0]);
c3 = curve2([-1 0 1],[0 -1/6 0]);
c  = geomcoerce('solid',{c1,c2,c3});

options = elOptionsSet('EPType','YModulPRatio','BCType', ...
   {'Dirichlet' 'Neumann'});

%Specify elasticity
%fn.lambda  = 3;
%fn.mu      = 1;
fn.YModul = 2;
fn.PRatio = 0.3;

%Specify the viscosity coefficient.
%fn.ViscCoef = {'viscCoef'};
fn.ViscCoef = 10;
fn.TimeStep = 1;

%Specify the boundary traction force.
fn.BndTFx = {[] 0};
fn.BndTFy = {[] 0};

%Specify boundary displacement.
fn.BndDispx = {0 []};
fn.BndDispy = {0 []};

%Myosin draging force.
fn.MyoDFx = {'myoDragFx' 0};
fn.MyoDFy = {'myoDragFy' 0};

fem.options = options;
fem.fn      = fn;
fem.fp      = [];

fem.geom    = c;
%fem.border  = 0;

%Group subdomains and boundaries.
fem.equ.ind = [1 2];
fem.bnd.ind = [1 2];

fem = elModelAssemble(fem);
fem = elasticSolve(fem,[]);
