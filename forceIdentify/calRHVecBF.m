%Calculate the right-hand side vector of the system for identifying the body
% force. 
% This vector is given by substracting from the data the solution to the 
% elastic equation with zero body force and the inhomogeous Dirichlet boundary  
% condition where the boundary displacement is provided by the data.
% It can only run after running 'setupModel' and loading data.

fprintf(1,'Constructing the right-hand vector : '); 

localStartTime = cputime;

%Set body force to be zero.
fp.BodyFx = {{'x' 'y'} {[] 0}};
fp.BodyFy = {{'x' 'y'} {[] 0}};

%Specify boundary condition.
for k = 1:numEdges
   fn.BndDispx{k} = 'bndDisp';
   fn.BndDispy{k} = 'bndDisp';
end
fem = elModelUpdate(fem,'fn',fn,'fp',fp);

rightU = cell(numTimeSteps,1);
for jj = 1:numTimeSteps
   for k = 1:numEdges
      fp.BndDispx{k} = {{'s'} {edgePPx{jj,k}}};
      fp.BndDispy{k} = {{'s'} {edgePPy{jj,k}}};
   end
   fem = elModelUpdate(fem,'fp',fp);

   fem = elasticSolve(fem,[]);

   [bspU1 bspU2] = postinterp(fem,'u1','u2',[dataPx{jj} dataPy{jj}].');

   rightU{jj} = zeros(2*numDP(jj),1);
   rightU{jj}(1:numDP(jj)) = dataU1{jj}-bspU1.';
   rightU{jj}(numDP(jj)+1:2*numDP(jj)) = dataU2{jj}-bspU2.';

   %For debugging, calculate the boundary displacements to compare with 
   % 'edge1U1' etc.
   for k = 1:numEdges
      [edgeUC1{jj,k} edgeUC2{jj,k}] = postinterp(fem,'u1','u2',edgeP{k});
   end
end

fprintf(1,'%f sec.\n', cputime-localStartTime);
