%Calculate the right-hand side vector of the system for the boundary traction
% force. 
% This vector is given by substracting from the data the solution to the 
% elastic equation with given body force (already identified), 
% free boundary condition on the edge where the boundary 
% traction force is to be identified and the inhomogeous Dirichlet boundary 
% condition on the other edges where the boundary displacement is provided by 
% the data.
% It can only run after running 'setupModel' and loading data and after the
% body force is identified.

fprintf(1,'\nConstructing the right-hand vector : ');

localStartTime = cputime;

%Set the body force.
load([resultPath 'bfId']);
fn.BodyFx = 'spMyoDFx';
fn.BodyFy = 'spMyoDFy';
fp.BodyFx = {{'x' 'y'} {fs coefBF(1:end/2)}};
fp.BodyFy = {{'x' 'y'} {fs coefBF(end/2+1:end)}};

for k = 1:numEdges
   BCTypes{k}     = 'Dirichlet';
   fn.BndDispx{k} = 'bndDisp';
   fn.BndDispy{k} = 'bndDisp';
   fp.BndDispx{k} = {{'s'} {edgePPx{k}}};
   fp.BndDispy{k} = {{'s'} {edgePPy{k}}};

   fn.BndTracFx{k} = 0;
   fn.BndTracFy{k} = 0;
end

rightU  = cell(numEdges,1);
for k = 1:numEdges
   BCTypes{k} = 'Neumann';
   options    = elOptionsSet(options,'BCType',BCTypes);

   fem = elModelUpdate(fem,'options',options,'fn',fn,'fp',fp);
   fem = elasticSolve(fem,[]);

   %Substract the solution from the data to get the right hand vector.
   [bspU1 bspU2] = postinterp(fem,'u1','u2',[dataPx dataPy].');

   rightU{k} = zeros(2*numDP,1);
   rightU{k}(1:numDP) = dataU1-bspU1;
   rightU{k}(numDP+1:2*numDP) = dataU2-bspU2;

   BCTypes{k} = 'Dirichlet';

   %For debugging, calculate the boundary displacements to compare with 'edgeU1'
   % and 'edgeU2'.
   for j = 1:numEdges
      [edgeUC1{j} edgeUC2{j}] = postinterp(fem,'u1','u2',edgeP{k});
   end

   fprintf(1,'   Edge %d (out of %d) finished.  %f sec.\n', ...
      k,numEdges,cputime-localStartTime);
end

fprintf(1,'Total time spent : %f sec.\n', cputime-localStartTime);
