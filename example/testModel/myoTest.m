clear all;

%Construct the geometry first.
c1 = curve2([-1 1],[-1 -1]);
c2 = curve2([-1 1],[1 1]);
c3 = curve2([-1 -1],[-1 1]);
c4 = curve2([1 1],[-1 1]);
c  = geomcoerce('solid',{c1,c2,c3,c4});

options = elOptionsSet('EPType','YModulPRatio','BCType', ...
   {'Dirichlet' 'Neumann'});

%Specify elasticity
%fn.lambda  = 3;
%fn.mu      = 1;
fn.YModul = 2;
fn.PRatio = 0.3;

%Specify the viscosity coefficient.
%fn.ViscCoef = {'viscCoef'};
fn.ViscCoef = 1;
fn.TimeStep = 1;

%Specify the boundary traction force.
fn.BndTFx = {[] 0};
fn.BndTFy = {[] 0};

%Specify boundary displacement.
fn.BndDispx = {0 []};
fn.BndDispy = {0 []};

%Myosin draging force.
fn.MyoDFx = {'testFx'};
fn.MyoDFy = {'testFy'};

fem.options = options;
fem.fn      = fn;
fem.fp      = [];

fem.geom    = c;
%fem.border  = 0;

%Group subdomains and boundaries.
fem.equ.ind = [1];
fem.bnd.ind = [1 2 2 1];

fem = elModelAssemble(fem);
fem = elasticSolve(fem,[]);
