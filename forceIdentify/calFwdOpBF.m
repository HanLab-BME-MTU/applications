%This script file computes the matrix approximation to the forward linear
% operator for the body force. It can only run after running 'setupModel' and 
% loading data.

fprintf(1,'\nConstructing the matrix A for the Body Force :\n');

localStartTime = cputime;

%Step 1: Build the basis of the function space for external body force using 
% the shape function in FEMLAB. We created a structure 'fs' to store 
% information about this function space. 'fs' is an FEM structure with 
% additional field.
%
%First, we need to set up the geometry and meshing of the domain for the
% external body force function space.
if isempty(bfDomPGx) 
   %If 'bfDomPG' is not defined, we recover the body force for the whole
   % region.
   bfDomPGx = recPGx;
   bfDomPGy = recPGy;
end

curvL = [bfDomPGx(bfDomPGVI(1):bfDomPGVI(2)) ...
         bfDomPGy(bfDomPGVI(1):bfDomPGVI(2))].';
curvT = [bfDomPGx(bfDomPGVI(2):bfDomPGVI(3)) ...
         bfDomPGy(bfDomPGVI(2):bfDomPGVI(3))].';
curvR = [bfDomPGx(bfDomPGVI(3):bfDomPGVI(4)) ...
         bfDomPGy(bfDomPGVI(3):bfDomPGVI(4))].';
curvB = [bfDomPGx(bfDomPGVI(4):end) bfDomPGy(bfDomPGVI(4):end)].';

%Name of any function in the space.
%fs.dim = {'v1' 'v2'};
fs.dim = {'v'};

fs.shape = 1;
fs.bnd.h = 1;
fs.bnd.r = 0;
fs.mesh  = rectmesh(curvL,curvT,curvR,curvB,bfDomHin,bfDomVin);
fs.geom  = fs.mesh;
fs.xmesh = meshextend(fs);

%The Degree of Freedom vector for the basis or shape function in the finite
% element space.
%dimBF  = flngdof(fs)/2; %Dimention of external Body Force.
dimBF  = flngdof(fs); %Dimention of external Body Force.
coefBF = zeros(dimBF,1);

%Step 2: Construct the matrix approximation to the forward operator.
%We use each basis function as the body force to solve our
% continuum mechanics system. The solution gives us each column of the matrix.
A   = zeros(numDP,2,dimBF,2);
sol = cell(2*dimBF,1); %To store the solutions corresponding to each basis.

col = 0;
ll  = 0;
if strcmp(fwdOpComputed,'none') == 1
   for k = 1:dimBF
      ll = ll+1;
      coefBF(k) = 1;
      fp.BodyFx = {{'x' 'y'} {fs coefBF}};
      fp.BodyFy = {{'x' 'y'} {[] 0}};
      fem = elModelUpdate(fem,'fp',fp);
      fem = elasticSolve(fem,[]);
      sol{col+1} = fem.sol;

      %'bspU1' and 'bspU2' is used to temporarily store the solution at the
      % data points.
      [bspU1 bspU2] = postinterp(fem,'u1','u2',[dataPx dataPy].');
      A(:,1,ll,1)  = bspU1.';
      A(:,2,ll,1)  = bspU2.';

      fp.BodyFx = {{'x' 'y'} {[] 0}};
      fp.BodyFy = {{'x' 'y'} {fs coefBF}};
      fem = elModelUpdate(fem,'fp',fp);
      fem = elasticSolve(fem,[]);
      sol{col+2} = fem.sol;

      [bspU1 bspU2] = postinterp(fem,'u1','u2',[dataPx dataPy].');
      A(:,1,ll,2)  = bspU1.';
      A(:,2,ll,2)  = bspU2.';

      coefBF(k) = 0;

      col = col+2;
      if rem(col,10) == 0
         fprintf(1,sprintf('   Columns %d out of %d finished.  %f sec.\n', ...
            col,2*dimBF,cputime-localStartTime));
      end
   end
   save([resultPath 'AbfId'],'A');
   save([modelPath 'solBFId'],'sol');
elseif strcmp(fwdOpComputed,'fem') == 1
   load([modelPath 'solBFId'],'sol');
   for k = 1:dimBF
      ll = ll+1;

      %'bspU1' and 'bspU2' is used to temporarily store the solution at the
      % data points.
      fem.sol = sol{col+1};
      [bspU1 bspU2] = postinterp(fem,'u1','u2',[dataPx dataPy].');
      A(:,1,ll,1)  = bspU1.';
      A(:,2,ll,1)  = bspU2.';

      fem.sol = sol{col+2};
      [bspU1 bspU2] = postinterp(fem,'u1','u2',[dataPx dataPy].');
      A(:,1,ll,2)  = bspU1.';
      A(:,2,ll,2)  = bspU2.';

      col = col+2;
      if rem(col,100) == 0
         fprintf(1,sprintf('   Columns %d out of %d finished.  %f sec.\n', ...
            col,2*dimBF,cputime-localStartTime));
      end
   end
   save([resultPath 'AbfId'],'A');
elseif strcmp(fwdOpComputed,'A') == 1
   load([resultPath 'AbfId'],'A');
end

fprintf(1,'   Total time spent : %f sec.\n', cputime-localStartTime);

