%Post assemble of the identified body force and some calculation
% for debugging purpose. It is run after 'solveLSBF'.

fprintf(1,'Post assemble of the identified body force : ');

localStartTime = cputime;

coefBFx = zeros(dimFS,numTimeSteps);
coefBFy = zeros(dimFS,numTimeSteps);
for jj = 1:numTimeSteps
   coefBFx(indDomDOF,jj) = coef(1:dimBF,jj);
   coefBFy(indDomDOF,jj) = coef(dimBF+1:2*dimBF,jj);
end
save([resultPath 'coefBFId'],'fs','coefBFx','coefBFy');

%For debugging : Solve the elastic equation with the identified force and
% compare the computed displacement with the measured displacement data.

%Specify boundary condition.
for k = 1:numEdges
   fn.BndDispx{k} = 'bndDisp';
   fn.BndDispy{k} = 'bndDisp';
end

fn.BodyFx = 'femBodyF';
fn.BodyFy = 'femBodyF';

fem = elModelUpdate(fem,'fn',fn,'fp',fp);

%Set the points where the identified force is to be calculated for
% demonstration.
if strcmp(forcePointType,'grid') == 1
   forcePx = gridPx{1};
   forcePy = gridPy{1};
else strcmp(forcePointType,'data') == 1
   forcePx = dataPx{1};
   forcePy = dataPy{1};
end

%Get the points that are inside the identification region.
[is,pe] = postinterp(fem,[forcePx forcePy].');
forcePx(pe) = [];
forcePy(pe) = [];

bodyFRx = cell(numTimeSteps,1);
bodyFRy = cell(numTimeSteps,1);
resDisp = cell(numTimeSteps,1);
resBF   = cell(numTimeSteps,1);
dataUC1 = cell(numTimeSteps,1);
dataUC2 = cell(numTimeSteps,1);
for jj = 1:numTimeSteps
   fp.BodyFx = {{'x' 'y'} {fs coefBFx(:,jj)}};
   fp.BodyFy = {{'x' 'y'} {fs coefBFy(:,jj)}};

   for k = 1:numEdges
      fp.BndDispx{k} = {{'s'} {edgePPx{jj,k}}};
      fp.BndDispy{k} = {{'s'} {edgePPy{jj,k}}};
   end

   fem = elModelUpdate(fem,'fp',fp);
   fem = elasticSolve(fem,[]);

   %Calculate the identified force.
   [bodyFRx{jj},bodyFRy{jj}] = postinterp(fem,'f1','f2', ...
      [forcePx forcePy].');
   %[bfxR,bfyR] = postinterp(fem,'f1','f2',fs.mesh.p);

   %The displacements computed with the identified force. To be compared with 
   % 'dataU1' and 'dataU2'.
   [dataUC1{jj} dataUC2{jj}] = postinterp(fem,'u1','u2', ...
      [dataPx{jj} dataPy{jj}].');

   %The displacements on the points where the force field is to be
   % demonstrated.
   [dispRU1{jj},dispRU2{jj}] = postinterp(fem,'u1','u2', ...
      [forcePx forcePy].');
   %[gridUC1 gridUC2] = postinterp(fem,'u1','u2',[gridPx gridPy].');

   %Calculate the residue of the forword computed displacements.
   resDisp{jj} = (dataUC1{jj}-dataU1{jj}).^2 + (dataUC2{jj}-dataU2{jj}).^2;

   %Calculate the residue of the regularized force.
   resBF{jj} = sigma*coef(:,jj).*coef(:,jj);

   %The displacement on the edge to be compared with 'edgeU1' etc.
   for k = 1:numEdges
      [edgeUC1{jj,k} edgeUC2{jj,k}] = postinterp(fem,'u1','u2',edgeP{k});
   end

   %Save the identified body force calculated on the demonstration points.
   save([resultPath 'bfId'],'forcePx','forcePy', ...
      'bodyFRx','bodyFRy','dispRU1','dispRU2');

   %Save the computed displacement from the identified body force.
   save([resultPath 'dispId'],'dataPx','dataPy','dataUC1','dataUC2');
end

fprintf(1,'%f sec.\n',cputime-localStartTime);
