clear options fn fp c fem;

options = elOptionsSet('EPType','YModulPRatio');

%fn.lambda  = 100;
%fn.mu      = 1;
fn.YModul = 2;
fn.PRatio = 0.2;
%fn.ViscCoef = {'viscCoef01'};
fn.ViscCoef = 0; %{'viscCoef01'};
fn.TimeStep = 1;
fn.BndDispx = {'circDispx'};
fn.BndDispy = {'circDispy'};

fp.BndDispy = {{1}};

c = circ2(0,0,1);

fem.options = options;
fem.fn      = fn;
fem.fp      = fp;

fem.geom    = c;
fem.equ.ind = [1];
fem.bnd.ind = [1 1 1 1];

fem = elModelAssemble(fem);
fem = elasticSolve(fem,[]);
