%Post assemble of the identified body force and some calculation
% for debugging purpose. It is run after 'solveLSBF'.

fprintf(1,'Post assemble of the identified body force : ');

localStartTime = cputime;

load([resultPath 'coefBFId']);

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
if strcmp(bfDisplaySite,'grid') == 1
   bfDisplayPx = gridPx{1};
   bfDisplayPy = gridPy{1};
elseif strcmp(bfDisplaySite,'data') == 1
   bfDisplayPx = dataPx{1};
   bfDisplayPy = dataPy{1};
end

%Get the points that are inside the identification region.
[is,pe] = postinterp(fem,[bfDisplayPx bfDisplayPy].');
bfDisplayPx(pe) = [];
bfDisplayPy(pe) = [];

recBF       = zeros(length(bfDisplayPx),2,numTimeSteps);
residueDisp = cell(numTimeSteps,1);
residueBF   = cell(numTimeSteps,1);
recDispU    = zeros(length(bfDisplayPx),2,numTimeSteps);
dataUC1     = cell(numTimeSteps,1);
dataUC2     = cell(numTimeSteps,1);
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
   [recBFx,recBFy] = postinterp(fem,'f1','f2', ...
      [bfDisplayPx bfDisplayPy].');
   recBF(:,:,jj) = [recBFx; recBFy].';
   %[bfxR,bfyR] = postinterp(fem,'f1','f2',fs.mesh.p);

   %The displacements on the points where the force field is to be
   % demonstrated.
   [recDispU1,recDispU2] = postinterp(fem,'u1','u2', ...
      [bfDisplayPx bfDisplayPy].');
   recDispU(:,:,jj) = [recDispU1;recDispU2].';

   %The displacements computed with the identified force. To be compared with 
   % 'dataU1' and 'dataU2'.
   [dataUC1{jj} dataUC2{jj}] = postinterp(fem,'u1','u2', ...
      [dataPx{jj} dataPy{jj}].');
   dataUC1{jj} = dataUC1{jj}.';
   dataUC2{jj} = dataUC2{jj}.';

   %Calculate the residue of the forword computed displacements.
   residueDisp{jj} = (dataUC1{jj}-dataU1{jj}).^2 + (dataUC2{jj}-dataU2{jj}).^2;

   %Calculate the residue of the regularized force.
   residueBF{jj} = sigma*coef(:,jj).*coef(:,jj);

   %The displacement on the edge to be compared with 'edgeU1' etc.
   for k = 1:numEdges
      [edgeUC1{jj,k} edgeUC2{jj,k}] = postinterp(fem,'u1','u2',edgeP{k});
   end

   %Save the identified body force calculated on the demonstration points.
   save([resultPath 'bfId'],'bfDisplayPx','bfDisplayPy', ...
      'recBF','recDispU','residueDisp','residueBF');

   %Save the computed displacement from the identified body force.
   save([resultPath 'dispId'],'dataPx','dataPy','dataUC1','dataUC2');
end

fprintf(1,'%f sec.\n',cputime-localStartTime);
