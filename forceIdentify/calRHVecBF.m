%Calculate the right-hand side vector of the system for identifying the body
% force. 
% This vector is given by substracting from the data the solution to the 
% elastic equation with zero body force and the inhomogeous Dirichlet boundary  
% condition where the boundary displacement is provided by the data.
% It can only run after running 'setupModel' and loading data.

fprintf(1,'Constructing the right-hand vector : '); 

localStartTime = cputime;

%Specify boundary condition.
for k = 1:numEdges
   fn.BndDispx{k} = 'bndDisp';
   fn.BndDispy{k} = 'bndDisp';
   fp.BndDispx{k} = {{'s'} {edgePPx{k}}};
   fp.BndDispy{k} = {{'s'} {edgePPy{k}}};
end

%Set body force to be zero.
fp.BodyFx = {{'x' 'y'} {[] 0}};
fp.BodyFy = {{'x' 'y'} {[] 0}};

fem = elModelUpdate(fem,'fn',fn,'fp',fp);
fem = elasticSolve(fem,[]);

[bspU1 bspU2] = postinterp(fem,'u1','u2',[dataPx dataPy].');

rightU = zeros(2*numDP,1);
rightU(1:numDP) = dataU1-bspU1;
rightU(numDP+1:2*numDP) = dataU2-bspU2;

%For debugging, calculate the boundary displacements to compare with 'edge1U1'
% etc.
for k = 1:numEdges
   [edgeUC1{k} edgeUC2{k}] = postinterp(fem,'u1','u2',edgeP{k});
end

fprintf(1,'%f sec.\n', cputime-localStartTime);
