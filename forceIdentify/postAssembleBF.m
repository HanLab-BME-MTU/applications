%Post assemble of the identified body force and some calculation
% for debugging purpose. It is run after 'solveLSBF'.

fprintf(1,'Post assemble of the identified body force : ');

localStartTime = cputime;

coefBFx = coef(1:dimBF);
coefBFy = coef(dimBF+1:2*dimBF);
save([resultPath 'bfId'],'fs','coefBFx','coefBFy');

%For debugging : Solve the elastic equation with the identified force and
% compare the computed displacement with the measured displacement data.

%Specify boundary condition.
for k = 1:numEdges
   fn.BndDispx{k} = 'bndDisp';
   fn.BndDispy{k} = 'bndDisp';
   fp.BndDispx{k} = {{'s'} {edgePPx{k}}};
   fp.BndDispy{k} = {{'s'} {edgePPy{k}}};
end

fn.BodyFx = 'femBodyF';
fn.BodyFy = 'femBodyF';
fp.BodyFx = {{'x' 'y'} {fs coefBFx}};
fp.BodyFy = {{'x' 'y'} {fs coefBFy}};

fem = elModelUpdate(fem,'fn',fn,'fp',fp);
fem = elasticSolve(fem,[]);

%Get the grid points that are inside the identification region.
[is,pe] = postinterp(fem,[gridPx gridPy].');
gridPx(pe) = [];
gridPy(pe) = [];

%Calculate the identified force on the grid points.
[gridBFxR,gridBFyR] = postinterp(fem,'f1','f2',[gridPx gridPy].');
[dataBFxR,dataBFyR] = postinterp(fem,'f1','f2',[dataPx dataPy].');

%Identified force:
[bfxR,bfyR] = postinterp(fem,'f1','f2',fs.mesh.p);

%The displacements computed with the identified force. To be compared with 
% 'dataU1' and 'dataU2'.
[dataUC1 dataUC2] = postinterp(fem,'u1','u2',[dataPx dataPy].');
[gridUC1 gridUC2] = postinterp(fem,'u1','u2',[gridPx gridPy].');

%Save the computed displacement from the identified body force.
save([resultPath 'dispId'],'dataPx','dataPy','dataUC1','dataUC2');

%The displacement on the edge to be compared with 'edgeU1' etc.
for k = 1:numEdges
   [edgeUC1{k} edgeUC2{k}] = postinterp(fem,'u1','u2',edgeP{k});
end

fprintf(1,'%f sec.\n',cputime-localStartTime);
