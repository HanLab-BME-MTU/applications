clear fem

%Create a unit circle centered at the origin with the command circ2:
fem.geom = circ2;

%Visualize the geometry by entering the command:
geomplot(fem)

%create a triangular mesh on the geometry defined in the fem.geom field:
fem.mesh = meshinit(fem);

%Visualize the resulting mesh:
meshplot(fem)

%function definition:
fem.functions{1}.type = 'inline';
fem.functions{1}.name = 'testfun(x,y)';
fem.functions{1}.expr = 'x.^2+y.^2';
%fem.functions{1}.dexpr = {'diff(sin(x)*cos(y),x)', ... 'diff(sin(x)*cos(y),y)'}

%Specify the PDE coefficients, using the coefficient form of the basic
%equations and set both c and f equal to 1 (see the section ?Using the
%Coefficient Form PDEs? on page 247 of the COMSOL Multiphysics Modeling
%Guide for details on the PDE terminology in COMSOL Multiphysics). All
%boundaries should have u = 0 as boundary conditions, which means setting h
%to 1. All coefficients you do not specify are zero by default. 
fem.equ.f = 'testfun(x,y)';
fem.equ.c = 1;
fem.bnd.h = 1;

%The default name of the dependent variable is u, but it is possible to use
%another name, for example T, by entering fem.dim = 'T'; however, in this
%example you do not have to specify fem.dim.

fem.dim = 'u';

%The PDE coefficients, coefficients in the boundary conditions, and initial
%conditions are not restricted to constant values. By specifying them as
%strings, you can define model parameters to be spatially varying (function
%of x, y, and z), time dependent (function of t), complex-valued (function
%of i or j), nonlinear (function of any and all dependent variables), and
%discontinuous. Furthermore, the string can reference an interpolation
%table or be a subroutine (any M-file in the current path or built-in
%functions in MATLAB). For more information about variable types that you
%can use and create, see ?Variables and Functions? on page 18. For example,
%you can define an initial condition that varies with x as
%atan(cos(0.5?x)):
%fem.init = 'atan(cos(0.5*pi*x))';
%fem.init = 'testfun(x,y)*x';


%Choose quadratic Lagrange elements:
fem.shape = 2;

%This setting must be made explicitly because the default Lagrange element
%order on the command line is 1.

%The meshextend function creates the extended mesh object, which is
%required for assembling the problem:
fem.xmesh = meshextend(fem);

%Solve the PDE and plot the solution:
fem.sol = femstatic(fem);
postplot(fem,'tridata','u','triz','u');

%Compute the maximal error by evaluating the difference with the exact
%solution:
pd = posteval(fem,'u-(1-x^2-y^2)/4');
er = max(max(pd.d))

%If the error is not sufficiently small, refine the mesh and update the
%model: 
fem.mesh = meshrefine(fem);
fem.xmesh = meshextend(fem);

%You can then solve the problem on the new mesh, plot the solution, and
%recompute the error by repeating Steps 7 and 8. 