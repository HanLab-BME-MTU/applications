clear options fn fp c fem;

options = elOptionsSet('EPType','YModulPRatio','BCType',{'Neumann'});

%fn.lambda  = 100;
%fn.mu      = 1;
fn.YModul = 2;
fn.PRatio = 0.2;
fn.ViscCoef = {'viscCoef'};
%fn.ViscCoef = 10;
fn.TimeStep = 1;
fn.BndTFx = 0;
fn.BndTFy = 0;

fn.MyoDFx = {'myoDragFx'};
fn.MyoDFy = {'myoDragFy'};

fp = [];

c = circ2(0,0,1);

fem.options = options;
fem.fn      = fn;
fem.fp      = fp;

fem.geom    = c;
fem.equ.ind = [1];
fem.bnd.ind = [1 1 1 1];

fem = elModelAssemble(fem);
fem = elasticSolve(fem,[]);
