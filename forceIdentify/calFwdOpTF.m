%This script file computes the matrix approximation to the forward linear
% operator for the body force. It can only run after running 'setupModel' and 
% loading data.
% For the function space of the boundary traction force, we use B-spline since
% it it one dimention.

%edgBrksTF : The breaks on each edge for the construction of B-spline.
%edgKnotsT : FThe knot sequence on each edge built from the breaks.
%dimTF     : Dimention of boundary Traction Force.
%bspTF     : The sp-form of the boundary Traction Force.
%coefTF    : The coefficient for the construction of the boundary Traction 
%            Force.
%A         : The matrix respresentation of the forward operator for the
%            boundary Traction Force.

fprintf(1,['\nConstructing the matrix A for the boundary ' ...
   'traction force :\n');

localStartTime = cputime;

edgBrksTF  = cell(numEdges,1);
edgKnotsTF = cell(numEdges,1);
dimTF      = zeros(numEdges,1); 
bspTF      = cell(numEdges,1);
coefTF     = cell(numEdges,1);
A          = cell(numEdges,1);
sol        = cell(numEdges,1);

for k = 1:numEdges
   numEdgBrksTF  = floor(edgArcLen(k)/edgBrkDistTF);
   edgBrksTF{k}  = linspace(0,edgArcLen(k),numEdgeBrksTF);
   edgKnotsTF{k} = augknt(edgBrksTF{k},bspOrderTF);
   dimTF(k)      = length(edgKnotsTF{k})-bspOrderTF;

   coefTF{k} = zeros(1,dimTF(k));
   bspTF{k}  = cell(dimTF(k),1);
end

if strcmp(fwdOpComputed,'none') == 1
   for k = 1:numEdges
      %Construct the sequence of breaks according to 'EdgBrkDistTF', the
      % distance between the breaks.
      A{k}   = zeros(numDP,2,dimTF{k},2);
      sol{k} = cell(2*dimTF(k),1);

      BCTypes{k} = 'Neumann';
      options    = elOptionsSet(options,'BCType',BCTypes);

      fem = elModelUpdate(fem,'options',options);
      for j = 1:dimTF(k)
         coefTF{k}(j) = 1;
         bspTF{k}{j} = spmak(knotS{k},coefTF{k});
         coefTF{k}(j) = 0;

         fp.BndTracFx{k} = {{'s'} {bspTF{k}{j}}};
         fp.BndTracFy{k} = {{'s'} {0}};
         fem = elModelUpdate(fem,'fp',fp);
         fem = elasticSolve(fem,[]);
         sol{k}{2*j-1} = fem.sol;

         %'bspU1' and 'bspU2' is used to temporarily store the solution at the
         % data points.
         [bspU1 bspU2] = postinterp(fem,'u1','u2',[dataPx dataPy].');
         A{k}(:,1,j,1)  = bspU1.';
         A{k}(:,2,j,1)  = bspU2.';

         fp.BndTracFx{k} = {{'s'} {0}};
         fp.BndTracFy{k} = {{'s'} {bspTF{k}{j}}};
         fem = elModelUpdate(fem,'fp',fp);
         fem = elasticSolve(fem,[]);
         sol{k}{2*j} = fem.sol;

         %'bspU1' and 'bspU2' is used to temporarily store the solution at the
         % data points.
         [bspU1 bspU2] = postinterp(fem,'u1','u2',[dataPx dataPy].');
         A{k}(:,1,j,2)  = bspU1.';
         A{k}(:,2,j,2)  = bspU2.';
      end
      fprintf(1,'   Edge %d (out of %d) finished.  %f sec.\n', ...
         k,numEdges,cputime-localStartTime);

      BCTypes{k} = 'Dirichlet';
   end

   save([resultPath 'AtfId'],'A');
   save([modelPath 'solTFId'],'sol');
elseif strcmp(fwdOpComputed,'fem') == 1
   load([modelPath 'solTFId'],'sol');
   for k = 1:numEdges
      %Construct the sequence of breaks according to 'EdgBrkDistTF', the
      % distance between the breaks.
      A{k}    = zeros(numDP,2,dimTF{k},2);

      for j = 1:dimTF(k)
         %'bspU1' and 'bspU2' is used to temporarily store the solution at the
         % data points.
         fem.sol = sol{k}{2*j-1};
         [bspU1 bspU2] = postinterp(fem,'u1','u2',[dataPx dataPy].');
         A{k}(:,1,j,1) = bspU1.';
         A{k}(:,2,j,1) = bspU2.';

         fem.sol = sol{k}{2*j};
         [bspU1 bspU2] = postinterp(fem,'u1','u2',[dataPx dataPy].');
         A{k}(:,1,j,2)  = bspU1.';
         A{k}(:,2,j,2)  = bspU2.';
      end
      fprintf(1,'   Edge %d (out of %d) finished.  %f sec.\n', ...
         k,numEdges,cputime-localStartTime);
   end

   save([resultPath 'AtfId'],'A');
elseif strcmp(fwdOpComputed,'A') == 1
   load([resultPath 'AtfId'],'A');
end


