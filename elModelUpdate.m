function fem = elModelUpdate(fem,varargin)
%ELMODELUPDATE Update the FEM structure that contains all the information
%              about the elastic equations. 
%
% SYNOPSIS : 
%    fem = elModelUpdate(fem,'name1',value1,'name2',value2,...)
%
% INPUT :
%    fem : The FEM structure of FEMLAB that is used to contain information of
%       the equation, domain and boundary condition etc.
%       See FEMLAB.
%    -------------------------------------------------------------------------
%    The rest input of property/value pairs are used to update or change the
%    specification of the model. Allowed properties are listed below:
%
%    geom : A FEMLAB geometry object that defines the geometry of the domain.
%       It will be directly passed to 'fem': 
%                fem.geom = geom
%       See FEMLAB.
%
%    options: A structure whose fields define some properties of the PDE 
%       system or provide options for solving the system.
%       Please see ELOPTIONSSET for an explaination of the fields and on how to 
%       create/alter this structure.
%
%    fn : A structure whose fields are used to define the material properties,
%       parameters or coefficient functions in the PDE or boundary conditions. 
%       See ELMODELASSEMBLE.
%
%    fp : A structure whose fields are cell arrays that pass parameters
%       to those functions defined in 'fn'. The names of the fields are the 
%       same as those of 'fn'. If one field in 'fn' is an evaluation string
%       and a numeric constant, the corresponding field in 'fp' can be left 
%       undefined or defined as an empty matrix, [].
%       See ELMODELASSEMBLE.
%
%    ind : The FEMLAB Index Vector that is used to group subdomains that 
%       share the same physical properties. It will be passed to 'fem':
%                fem.ind = ind;
%       See FEMLAB.
%
%    bndInd : The FEMLAB Index Vector that is used to group boundaries that 
%       share the same physical properties. It will be passed to 'fem':
%                fem.bnd.ind = bndInd;
%       See FEMLAB.
%
%    mesh : Modify or redefine the meshing of the domain. The value is passed
%       to various meshing tool provided by FEMLAB.
%    --------------------------------------------------------------------------
%
% OUTPUT :
%    fem : The FEM structure in FEMLAB that upon a successful exit contain the 
%       updated information about the PDE model and the meshing.
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

%First, we check if the input is of legitimate data type and if the property
% names are valid.
if ~isstruct(fem)
   error(['The first argument must be a structure. More precisely, ' ...
      'the FEM structure of FEMLAB.']);
end

if ceil((nargin-1)/2) ~= floor((nargin-1)/2)
   error('The property/value inputs must appear in pair.');
end

%Initialize 'options', 'fn', 'fp', 'ind' and 'bndInd' with the values in the
% old 'fem'.
if isfield(fem,'geom')
   geom = fem.geom;
else
   geom = [];
end

if isfield(fem,'mesh')
   mesh = fem.mesh;
else
   mesh = [];
end

if isfield(fem,'options')
   options = fem.options;
else
   options = [];
end

if isfield(fem,'fn')
   fn = fem.fn;
else
   fn = [];
end

if isfield(fem,'fp')
   fp = fem.fp;
else
   fp = [];
end

if isfield(fem,'equ') & isfield(fem.equ,'ind')
   ind = fem.equ.ind;
else
   ind = [];
end

if isfield(fem,'bnd') & isfield(fem.bnd,'ind')
   bndInd = fem.bnd.ind;
else
   ind = [];
end

%Allowed propery names.
propertyNames = {'geom'
                 'options'
                 'fn'
                 'fp'
                 'ind'
                 'bndInd'
                 'mesh'};

ReAssemble   = 'no';
geomChanged  = 'no';
meshChanged  = 'no';
fnChanged    = 'no';
fpChanged    = 'no';
for k = 1:2:nargin-1
   if ~ischar(varargin{k})
      msg = [sprintf('Property/value pair no. %d does not start ', (k+1)/2) ...
         'with a string that specify the name of the property.'];
      error(msg);
   end
   switch varargin{k}
      case 'geom'
         geom = varargin{k+1};
         geomChanged = 'yes';
      case 'mesh'
         mesh = varargin{k+1};
         meshChanged = 'yes';
      case 'options'
         options = varargin{k+1};
         ReAssemble = 'yes';
      case 'fn'
         fn = varargin{k+1};
         fnChanged = 'yes';
      case 'fp'
         fp = varargin{k+1};
         fpChanged = 'yes';
      case 'ind'
         ind = varargin{k+1};
         ReAssemble = 'yes';
      case 'bndInd'
         bndInd = varargin{k+1};
         ReAssemble = 'yes';
      otherwise
         error(['''' varargin{k} ''' is not a valid property name.']);
   end
end

if strcmp(geomChanged,'yes') | strcmp(meshChanged,'yes')
   if strcmp(meshChanged,'yes')
      fem = elModelAssemble(geom,mesh,options,fn,fp,ind,bndInd);
   else
      fem = elModelAssemble(geom,[],options,fn,fp,ind,bndInd);
   end
   return;
end

if strcmp(meshChanged,'yes')
   if isstruct(mesh)
      %'mesh' is a completely defined MESH structure. We define 'geom' using 
      % mesh as geometry in this case.
      fem = elModelAssemble([],mesh,options,fn,fp,ind,bndInd);
   else
      %'mesh' is a cell array of meshing parameters.
      fem = elModelAssemble(geom,mesh,options,fn,fp,ind,bndInd);
   end
   return;
end

if strcmp(ReAssemble,'yes')
   fem = elModelAssemble(geom,mesh,options,fn,fp,ind,bndInd);
   return;
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
if strcmp(fnChanged,'yes') == 1
   fem = inputsCheck(fem,options,fn,fp,ind,bndInd);
   fem = varDefine(fem);
   fem.const = constDefine(fem);
   fem.xmesh = meshextend(fem);
elseif strcmp(fpChanged,'yes') == 1
   fem = inputsCheck(fem,options,fn,fp,ind,bndInd);
   fem.const = constDefine(fem);
end
