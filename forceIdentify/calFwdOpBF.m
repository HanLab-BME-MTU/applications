%This script file computes the matrix approximation to the forward linear
% operator for the body force. It can only run after running 'setupModel' and 
% 'calFwdOpBF' or loading data.

fprintf(1,'Calculating the forward operator for the Body Force :\n');

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
%Find the indices of the DOFs whose shape functions have support disconnected 
% from the boundary.
N = assemble(fs,'out','N');
indDomDOF = find(full(sum(N,1))==0);

%The Degree of Freedom vector for the basis or shape function in the finite
% element space.
%Dimention of the function space where the basis can be nonzero on the 
% boundary.
dimFS  = flngdof(fs); 
coefFS = zeros(dimFS,1);
%Dimention of the function space where the basis function is kept zero on the
% boundary.
dimBF = length(indDomDOF);

%Step 2: Construct the matrix approximation to the forward operator.
%We use each basis function as the body force to solve our
% continuum mechanics system. The solution gives us each column of the matrix.
col = 0;
if strcmp(fwdOpComputed,'none') == 1
   fprintf(1,'\n  Solving the elastic equation for each basis function :\n');
   %To store the solutions corresponding to each basis.
   sol = cell(2*dimBF,1); 

   for j = 1:dimBF
      %'k' is the index of 'coefFS' whose corresponding basis function is zero
      % on the boundary.
      k = indDomDOF(j);
      coefFS(k) = 1;
      fp.BodyFx = {{'x' 'y'} {fs coefFS}};
      fp.BodyFy = {{'x' 'y'} {[] 0}};
      fem = elModelUpdate(fem,'fp',fp);
      fem = elasticSolve(fem,[]);
      sol{col+1} = fem.sol;

      fp.BodyFx = {{'x' 'y'} {[] 0}};
      fp.BodyFy = {{'x' 'y'} {fs coefFS}};
      fem = elModelUpdate(fem,'fp',fp);
      fem = elasticSolve(fem,[]);
      sol{col+2} = fem.sol;

      coefFS(k) = 0;

      col = col+2;
      if rem(col,10) == 0
         fprintf(1,'    %d basis functions out of ', col);
         fprintf(1,'%d finished.  %f sec.\n', ...
            2*dimBF,cputime-localStartTime);
      end
   end
   save([modelPath 'solBFId'],'sol');
elseif strcmp(fwdOpComputed,'fem') == 1
   fprintf(1,['  Loading the solution to the elastic equation ' ...
      'for each basis function.\n']);
   load([modelPath 'solBFId'],'sol');
end

if strcmp(fwdOpComputed,'A') == 1
   fprintf(1,'  Loading the matrix A for the Body Force.\n');
   load([resultPath 'AbfId'],'A');
else
   A   = cell(numTimeSteps,1);
   for jj = 1:numTimeSteps
      fprintf(1,'  Constructing the matrix A at time step %d :\n', jj);
      A{jj} = zeros(numDP(jj),2,dimBF,2);

      col = 0;
      ll  = 0;
      for k = 1:dimBF
         ll = ll+1;

         %'bspU1' and 'bspU2' is used to temporarily store the solution at the
         % data points.
         fem.sol = sol{col+1};
         [bspU1 bspU2] = postinterp(fem,'u1','u2',[dataPx{jj} dataPy{jj}].');
         A{jj}(:,1,ll,1)  = bspU1.';
         A{jj}(:,2,ll,1)  = bspU2.';

         fem.sol = sol{col+2};
         [bspU1 bspU2] = postinterp(fem,'u1','u2',[dataPx{jj} dataPy{jj}].');
         A{jj}(:,1,ll,2)  = bspU1.';
         A{jj}(:,2,ll,2)  = bspU2.';

         col = col+2;
         if rem(col,100) == 0
            fprintf(1,'    Columns %d out of %d finished.  %f sec.\n', ...
               col,2*dimBF,cputime-localStartTime);
         end
      end
      %Reshape 'A';
      A{jj} = reshape(A{jj},2*numDP(jj),2*dimBF);
   end
   save([resultPath 'AbfId'],'A');
end

fprintf(1,'  Total time spent : %f sec.\n', cputime-localStartTime);

