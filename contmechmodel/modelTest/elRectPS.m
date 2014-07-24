%This script file can be used as an example that shows how we build the
% foward elastic model by passing parameters to 'elModelAssemble' and solve it.
% It is also used to test if 'elModelAssemble' is correctly programmed to set
% up the model by comparing the results with those calculated through 
% the Structural Mechanics Module provided by Femlab. The corresponding model
% file is 'femRectPS.m' which can be loaded through Femlab GUI.

%Build the geometry. It is a rectangular with a circle in the middle. The body
% force is confind inside the circle.
R1 = rect2(-200,200,-150,150);
E1 = circ2(0,0,50,0);
G1 = geomcoerce('solid',{R1,E1});

%Pass the geometry to the 'fem' structure.
fem.geom = G1;

%Group the boundary
bndInd = { [1 4] [2 3]};

%Define boundary conditions.
% Neumann condition (force load) on group 1 boundaries.
% 'nx' 'ny': femlab variables that gives the unit normal direction pointing
% outwards from the domain.
fn.BndTracFx = { '1e-2*nx' [] };
fn.BndTracFy = { '1e-2*ny' [] };

%Dirichilet condition (constraints) on group 2 boundaries.
fn.BndDispx = { [] '-2*nx' };
fn.BndDispy = { [] '-2*ny' };

%Group the subdomains.
ind = [1 2];

%Define the body force. Restricted to domain 2 in this case.
fn.BodyFx = {0 '1e-5*sqrt(x^2+y^2)'};
fn.BodyFy = {0 '-1e-5*sqrt(x^2+y^2)'};

%Set the elastic parameters
A = 2;
V = 1;
fn.YModul = 'testYMod';
fp.YModul = { {'x' 'y'} {A V} };
fn.PRatio = 0.33;

%Set the options.
options = elOptionsSet('BCType',{'Neumann' 'Dirichlet'},'EPType', ...
   'YModulPRatio');

%Pass parameters to 'elModelAssemble'.
fem = elModelAssemble(G1,[],options,fn,fp,ind,bndInd);

%Solve the elastic equation.
fem = elasticSolve(fem,[]);

%Plot the solution.
figure(gcf); clf;
geomplot(fem,'edgelabel','on','sublabel','on'); hold on;
postplot(fem,'arrowdata',{'u1' 'u2'},'arrowcolor','k','arrowscale',0.5, ...
   'tridata','sqrt(u1^2+u2^2)','tribar','on');
axis equal; hold off;
