function fem = elModelAssemble(geom,mesh,options,fn,fp,ind,bndInd)
%ELMODELASSEMBLE  Assemble information of the elastic equations into the FEM
%                 structure of FEMLAB. 
%
% SYNOPSIS : 
%    fem = elModelAssemble(geom,mesh,options,fn,fp,ind,bndInd)
%
% INPUT :
%    geom : A FEMLAB geometry object that defines the geometry of the domain.
%       It will be directly passed to 'fem'
%                fem.geom = geom;
%       if it is not empty. Otherwise, we will define the geometry with the
%       mesh as geometry concept in FEMLAB, i.e.
%                fem.geom = mesh;
%    mesh: A cell array of meshing parameters that can be passed directly to
%       the 'meshinit' function in FEMLAB, or a completely defined MESH 
%       structure in FEMLAB. It can also be the empty matrix, []. In this
%       case, however, 'geom' can not be empty and a default meshing will be 
%       used on the input 'geom'.
%       See FEMLAB.
%    -------------------------------------------------------------------------
%    options: A structure whose fields define some properties of the PDE 
%       system or provide options for solving the system.
%       Please see ELOPTIONSSET for an explaination of the fields and on how to 
%       create/alter this structure.
%    --------------------------------------------------------------------------
%    fn : A structure whose fields are used to define the material properties,
%       parameters or coefficient functions in the PDE or boundary conditions. 
%
%       FIELD     : INTERPRETATION
%       ---------------------------------------------------------------
%       lambda    : the Lame parameter, lambda.
%       mu        : the Lame parameter, mu.
%       YModul    : Young's Modulus.
%       PRatio    : Poisson's Ratio.
%       VDragCoef : Viscosity Dragging Coeffecient.
%       TimeStep  : The period of time over which the displacement of the actin
%                   meshwork is measured.
%       ActVx     : The x-component and
%       ActVy     : The y-component of the Actin meshwork flow Velocity.
%       MyoDragFx : The x-component and
%       MyoDragFy : The y-component of the Myosin Dragging Force.
%       BodyFx    : The x-component and
%       BodyFy    : The y-component of any other unknown external Body Force.
%       ActPolyR  : Actin Polymerization Rate.
%       ActDPolyR : Actin De-Polymerization Rate.
%       BndTracFx : The x-component and
%       BndTracFy : the y-component of the Boundary Traction Force or 
%                   external loading force on the leading edge. 
%       BndDispx  : The x-component and
%       BndDispy  : The y-component of the Boundary Displacement. 
%       Init1     : First and 
%       Init2     : Second Initial conditions for dynamic problem.
%                   If there is only one time derivative, only 'Init1' is used.
%                   If there is no time derivative, they are ignored.
%       ---------------------------------------------------------------
%       Each field is a cell array of evaluation strings, function names, 
%       function handles or numerical values. They can also be numerical 
%       vectors. The length of the cell array or numerical vector equals the 
%       number of subdomains in the case of a domain dependent variable and
%       the number of boundaries in the case of a boundary dependent variable. 
%
%       The number of subdomains or boundaries does not refer to the
%       number of physical subdomains (or boundaries) in the geometry
%       but refer to the number of groups of subdomains that share the
%       same physical properties. Such a grouping is done through the
%       Index Vector concept of FEMLAB. See 'fem.ind' and 'fem.bnd.ind'.
%
%       If there is only one subdomain or if you want the same function or 
%       value to be applied to all the subdomains, you can either use a cell 
%       array of length 1 or you can just specify the string (or the 
%       function handle or the value) without using a cell array. 
%    
%       The parameters to the function is provided by the structure 'fp'
%       explained below. The distinction between an evaluation string and a
%       function name is determined by checking if there is any arguments
%       passed through the structure 'fp' to the function. 
%       See explaination for 'fp'.
%    
%       The field can be left undefined or defined as an empty matrix, 
%       [], if it is not used. For evaluation string, you can only use
%       variables that can be recogonized by FEMLAB.
%
%       The name of the function can not be the same as the correponding field
%       name. For example,
%          fn.YModul = 'YModul'; %The name of the function that defines the
%                                % Young's Modulus.
%          fp.YModul = { {'x' 'y'} {2 4} }; %Arguments to be passed to the
%                                           % function.
%       will give you an error.
%
%       See examples.
%
%    --------------------------------------------------------------------------
%    fp : A structure whose fields are cell arrays that pass parameters
%       to those functions defined in 'fn'. The names of the fields are the 
%       same as those of 'fn'. If one field in 'fn' is an evaluation string
%       and a numeric constant, the corresponding field in 'fp' can be left 
%       undefined or defined as an empty matrix, [].
%
%       The parameters to be passed are always seperated into two groups given
%       by two cell arrays. The first group is always a cell array of
%       variables recognizable by FEMLAB such as 't', 'x', 'y', 's' etc. and 
%       any user defined variables. The second group is basically a cell array 
%       of the rest of the parameters. The two cell arrays are then put
%       together in one upper level cell array. If you have more than one
%       subdomains or boundaries, you will need one more level of cell array.
%
%       See examples.
%    -------------------------------------------------------------------------
%    ind : The FEMLAB Index Vector that is used to group subdomains that 
%       share the same physical properties. It will be passed to 'fem':
%                fem.ind = ind;
%       See FEMLAB.
%    -------------------------------------------------------------------------
%    bndInd : The FEMLAB Index Vector that is used to group boundaries that 
%       share the same physical properties. It will be passed to 'fem':
%                fem.bnd.ind = bndInd;
%       See FEMLAB.
%    --------------------------------------------------------------------------
%
% OUTPUT :
%    fem : The FEM structure in FEMLAB that is used to contain the information 
%       about the PDE model, the meshing and the solution (after calling
%       ELASTICSOLVE). On exit, it will have all the information that is 
%       needed by FEMLAB to solve the model.
%    -------------------------------------------------------------------------
%
% EXAMPLES of 'fn' and 'fp' :
%    1. Evaluation string and numerical values:
%       fn.lambda = 'x^2+1';
%       fn.lambda = 2;
%       fn.lambda = {'x*y+1' 2}; %There are two subdomains.
%
%       Note: Please leave 'fp.lambda' undefined or define it to be the
%       empty matrix, [].
%
%    2. Function name or handle.
%       %One subdomain or it is applied to all subdomains.
%       fn.mu = 'f_mu';
%       A_mu = 9;
%       fp.mu = {{'x' 'y'} {4 5 A_mu}}; 
%
%       The syntax for function 'f_mu' is then
%                y = f_mu(x,y,p1,p2,p3);
%       Note: 'x' and 'y' are FEMLAB reserved variables for space coordinates.
%
%       %Multiple subdomains.
%       fn.mu = {'f_mu1' @f_mu2 'x+y'}; %3 subdomains.
%       fp.mu = { {{'x' 'y'} {3}}, {[] {3 4}}, []}
%
%       The syntax for function 'f_mu1' and 'f_mu2' are then
%                y = f_mu1(x,y,p1); %and
%                y = f_mu2(p1,p2);
%       Since for the 3rd subdomain we have an evaluation string, 'x+y', use
%       the empty matrix, [], to distinguish it from a function name.

if isempty(geom)
   if isempty(mesh) | ~isa(mesh,'struct')
      error(['When the first argument,''geom'' is empty, ' ...
         'the second argument, ''mesh'' must be a complete predefined ' ...
         'MESH structure of FEMLAB.']);
   else
      fem.geom = mesh;
   end
else
   if ~isempty(mesh) & ~isa(mesh,'struct') & ~iscell(mesh)
      error(['The second argument ''mesh'' has to be either the ' ...
         'empty matrix, [], or a cell array of meshing arguments, ' ...
         'or the completely defined MESH structure.']); 
   end
   fem.geom = geom; % set geometry.
end

%INPUTSCHECK is a private function that checks the legitimacy of the inputs
% and set default values for fields of 'options' that are not defined by the 
% user. On a successful return, all the inputs will be bundled to the 
% structure, 'fem' in the following way:
%    fem.options = options;
%    fem.fn      = fn;
%    fem.fp      = fp;
%    fem.ind     = ind;
%    fem.bnd.ind = bndInd;
fem = inputsCheck(fem,options,fn,fp,ind,bndInd);

%Start the construction of the fields of FEM structure that are needed by
% FEMLAB.
fem.sdim = {'x','y'}; % independent variables
fem.dim  = {'u1','u2'}; % dependent or solution variables
fem.form = 'coefficient';

%Since each parameter given in the FN structure (e.g. 'lambda') have different
% definition on different subdomains (or boundaries) and FEMLAB does not
% support variables that represent vectors or cell arrays, we need to define
% one variable for each subdomain (or boundary). The name of these variables
% would be a combination of the parameter name (e.g. 'lambda') with the
% subdomain or boundary number (e.g. 'lambda1'). We designed a private
% function to accomplish this task.
fem.varNames = varNameDefine(fem);

%Define variables that can be used by FEMLAB to get all the coefficients and 
% forces in the PDE and all the functions that are needed in the boundary 
% conditions. 
fem = varDefine(fem);

%Define constants that can be used by FEMLAB to get constant coefficients in 
% the PDE.
fem = constDefine(fem);

%Get the name list for each variable from 'fem.varNames' for easy access.
% Therefore, instead of getting the name of the variable from
% 'fem.varNames.YModul{j}', we can use 'YModul{j}' where 'j' refers to the jth
% subdomain or boundary.
YModul    = fem.varNames.YModul; % e.g. YModul = {YModul1 YModul2}.
PRatio    = fem.varNames.PRatio;
lambda    = fem.varNames.lambda;
mu        = fem.varNames.mu;
VDragCoef = fem.varNames.VDragCoef; 
TimeStep  = fem.varNames.TimeStep; 
ActVx     = fem.varNames.ActVx;
ActVy     = fem.varNames.ActVy;
MyoDragFx = fem.varNames.MyoDragFx;
MyoDragFy = fem.varNames.MyoDragFy;
BodyFx    = fem.varNames.BodyFx;
BodyFy    = fem.varNames.BodyFy;
ActPolyR  = fem.varNames.ActPolyR;
ActDPolyR = fem.varNames.ActDPolyR;
BndTracFx = fem.varNames.BndTracFx;
BndTracFy = fem.varNames.BndTracFy;
BndDispx  = fem.varNames.BndDispx;
BndDispy  = fem.varNames.BndDispy;
Init1     = fem.varNames.Init1;
Init2     = fem.varNames.Init2;

%Define the coefficients and forces that are used by FEMLAB. Please refer to
% FEMLAB for the meaning of all the coefficients.
fem.equ.c  = cell(1,fem.numSubDoms);
fem.equ.f  = cell(1,fem.numSubDoms);
fem.equ.a  = cell(1,fem.numSubDoms);
fem.equ.al = cell(1,fem.numSubDoms);
fem.equ.ga = cell(1,fem.numSubDoms);
fem.equ.be = cell(1,fem.numSubDoms);
fx         = cell(1,fem.numSubDoms);
fy         = cell(1,fem.numSubDoms);
for k = 1:fem.numSubDoms
   %Set the coefficient, 'c'.
   fem.equ.c{k} = { {[lambda{k} '+2*' mu{k}] 0; 0 mu{k}} ...
                    {0 lambda{k}; mu{k} 0}; ...
                    {0 mu{k}; lambda{k} 0} ...
                    {mu{k} 0; 0 [lambda{k} '+2*' mu{k}]} };

   %Set the coefficient, 'a' which is determined by the viscosity coefficient,
   % 'fn.VDragCoef' and the time step, 'fn.TimeStep' over which the actin
   % meshwork has moved the displacement 'ux' and 'uy'. 
   if ~isempty(VDragCoef{k})
      if isempty(TimeStep{k})
         error(['''fn.TimeStep'' has to be defined when you have '...
            '''fn.VDragCoef'' defined.']);
      end

      fem.equ.a{k} = { [VDragCoef{k} './' TimeStep{k}] 0; ...
                       0 [VDragCoef{k} './' TimeStep{k}] };
   else
      fem.equ.a{k} = { 0 0; 0 0 };
   end

   %By FEMLAB default, all the other coefficients, 'al', 'ga' and 'be' might
   % be zero. But to make sure, we set them to be zero explicitly here.
   fem.equ.al{k} = { {0; 0} {0; 0}; {0; 0} {0; 0} };
   fem.equ.ga{k} = { {0; 0}; {0; 0} };
   fem.equ.be{k} = { {0; 0} {0; 0}; {0; 0} {0; 0} };

   %Set the force term, 'f'.
   fx{k} = '0';
   fy{k} = '0';
   if ~isempty(MyoDragFx{k}) 
      %isfield(fem.fn,'MyoDragFx') & ~isempty(fem.fn.MyoDragFx{k})
      fx{k} = [fx{k} '+' MyoDragFx{k}];
   end

   if ~isempty(MyoDragFy{k}) 
      %isfield(fem.fn,'MyoDragFy') & ~isempty(fem.fn.MyoDragFy{k})
      fy{k} = [fy{k} '+' MyoDragFy{k}];
   end
   
   if ~isempty(BodyFx{k}) 
      fx{k} = [fx{k} '+' BodyFx{k}];
   end

   if ~isempty(BodyFy{k}) 
      fy{k} = [fy{k} '+' BodyFy{k}];
   end

    fem.equ.f{k} = {fx{k}; fy{k}};
end


%Define the boundary condition
fem.bnd.h = cell(1,fem.numBnds);
fem.bnd.r = cell(1,fem.numBnds);
fem.bnd.g = cell(1,fem.numBnds);
fem.bnd.q = cell(1,fem.numBnds);
for k = 1:fem.numBnds
   if strcmp(fem.options.BCType{k},'Dirichlet') == 1
      fem.bnd.h{k} = {1 0; 0 1};
      fem.bnd.r{k} = {BndDispx{k}; BndDispy{k}};

      %Set 'q' to be zero although it might be zero by FEMLAB default.
      %fem.bnd.g{k} = {0; 0};
      %fem.bnd.q{k} = {0 0; 0 0};
   elseif strcmp(fem.options.BCType{k},'Neumann') == 1
      fem.bnd.h{k} = {0 0; 0 0};
      fem.bnd.r{k} = {0; 0};
      fem.bnd.g{k} = {BndTracFx{k}; BndTracFy{k}};

      %Set 'q' to be zero although it might be zero by FEMLAB default.
      fem.bnd.q{k} = {0 0; 0 0};
   end
end

%Define the shape function.
fem.shape = 2; % Largrange elelments.

%Define the meshing.
if isempty(mesh)
   fem.mesh = meshinit(fem);
   fem.mesh = meshrefine(fem);
   fem.mesh = meshrefine(fem);
elseif isa(mesh,'struct')
   fem.mesh = mesh;
   fem.geom = mesh;
else 
   %'mesh' is a cell array of meshing parameters that are accepted by
   % 'meshinit'.
   fem.mesh = meshinit(fem,mesh{:});
end

%Extending the mesh.
fem.xmesh = meshextend(fem);

%%%%%% End of Main Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

