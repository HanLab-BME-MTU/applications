%Post assemble of the identified body force and some calculation
% for debugging purpose. It is run after 'solveAdhLConstrnBF'.

%Get one cell image
cellImg = imread(imgFile{1});

if strcmp(fwdOpComputed,'all') == 1
   fprintf(1,'Load the identified body force and displacements: \n');
   load([resultPath 'bfAdhId']);
   load([resultPath 'dispAdhId']);
   return;
end

fprintf(1,'Post assemble of the identified contraction force and adhesion: ');

localStartTime = cputime;

load([modelPath 'femId']);
load([resultPath 'dispField']);
load([resultPath 'dataDisp']);
load([resultPath 'edgeDisp']);
load([resultPath 'coefAdhLConstrnBF']);

%For debugging : Solve the elastic equation with the identified force and
% compare the computed displacement with the measured displacement data.

%Specify boundary condition.
for k = 1:numEdges
   fn.BndDispx{k} = 'bndDisp';
   fn.BndDispy{k} = 'bndDisp';
end

fn.BodyFx = 'femBodyF';
fn.BodyFy = 'femBodyF';

fn.VDragCoef = 'femBodyF';
fn.TimeStep  = 1;

fem = elModelUpdate(fem,'fn',fn);

%Set the points where the identified force is to be calculated for
% demonstration.
if strcmp(bfDisplaySite,'grid') == 1
   bfDisplayPx = gridPx;
   bfDisplayPy = gridPy;
elseif strcmp(bfDisplaySite,'speckle') == 1
   bfDisplayPy = speckleP{1}(:,1);
   bfDisplayPx = speckleP{1}(:,2);
end

%Get the points that are inside the identification region.
[is,pe] = postinterp(fem,[bfDisplayPx bfDisplayPy].');
bfDisplayPx(pe) = [];
bfDisplayPy(pe) = [];

%Get the grid points that are inside the identification region. They will be
% used to calculate the intensity score for contraction and adhesion.
[is,pe] = postinterp(fem,[gridX gridY].');
gridIn = 1:length(gridX);
gridIn(pe) = [];

%The data points must also be inside the target recovery region.
for jj = 1:numTimeSteps
   [is,pe] = postinterp(fem,[dataPx{jj} dataPy{jj}].');
   dataPx{jj}(pe) = []; % 'pe': index of points outside 'msh'.
   dataPy{jj}(pe) = [];
   dataU1{jj}(pe) = [];
   dataU2{jj}(pe) = [];
end

recAdhMap   = zeros(size(cellImg,1),size(cellImg,2),numTimeSteps);
recMCF      = zeros(length(bfDisplayPx),2,numTimeSteps);
recADF      = zeros(length(bfDisplayPx),2,numTimeSteps);
residueDisp = cell(numTimeSteps,1);
residueBF   = cell(numTimeSteps,1);
recDispU    = zeros(length(bfDisplayPx),2,numTimeSteps);
dataUC1     = cell(numTimeSteps,1);
dataUC2     = cell(numTimeSteps,1);
dataMCFx    = cell(numTimeSteps,1);
dataMCFy    = cell(numTimeSteps,1);
dataAdh     = cell(numTimeSteps,1);
for jj = 1:numTimeSteps
   fp.BodyFx = {{'x' 'y'} {fs coefMCFx(:,jj)}};
   fp.BodyFy = {{'x' 'y'} {fs coefMCFy(:,jj)}};

   fp.VDragCoef = {{'x' 'y'} {fs coefAdh(:,jj)}};

   for k = 1:numEdges
      fp.BndDispx{k} = {{'s'} {edgePPx{jj,k}}};
      fp.BndDispy{k} = {{'s'} {edgePPy{jj,k}}};
   end

   fem = elModelUpdate(fem,'fp',fp);
   fem = elasticSolve(fem,[]);

   %Calculate the identified force.
   [recMCFx,recMCFy] = postinterp(fem,'f1','f2', ...
      [bfDisplayPx bfDisplayPy].');
   recMCF(:,:,jj) = [recMCFx; recMCFy].';
   %[bfxR,bfyR] = postinterp(fem,'f1','f2',fs.mesh.p);

   %The displacements on the points where the force field is to be
   % demonstrated.
   [recDispU1,recDispU2] = postinterp(fem,'u1','u2', ...
      [bfDisplayPx bfDisplayPy].');
   recDispU(:,:,jj) = [recDispU1;recDispU2].';

   %Calulate the adhesion dragging force on the demonstration points.
   recAdh   = postinterp(fem,'VDragCoef',[bfDisplayPx bfDisplayPy].');
   recADF(:,:,jj) = -[recAdh.*recDispU1;recAdh.*recDispU2].';

   %Produce the identified adhesion intensity map.
   [pixelH pixelV] = meshgrid([1:size(cellImg,2)],[1:size(cellImg,1)]);
   pixelx = reshape(pixelH,1,length(pixelH(:)));
   pixely = reshape(pixelV,1,length(pixelV(:)));
   [adhI,pe] = postinterp(fem,'VDragCoef',[pixelx;pixely]);
   adhI(pe)  = 0;
   recAdhMap(:,:,jj) = reshape(adhI,size(cellImg,1),size(cellImg,2));

   %The displacements computed with the identified force and the identified
   % force at the data points. To be compared with 
   % 'dataU1' and 'dataU2'.
   [dataUC1{jj} dataUC2{jj}] = postinterp(fem,'u1','u2', ...
      [dataPx{jj} dataPy{jj}].');
   dataUC1{jj} = dataUC1{jj}.';
   dataUC2{jj} = dataUC2{jj}.';

   [dataMCFx{jj} dataMCFy{jj}] = postinterp(fem,'f1','f2', ...
      [dataPx{jj} dataPy{jj}].');
   dataMCFx{jj} = dataMCFx{jj}.';
   dataMCFy{jj} = dataMCFy{jj}.';

   dataAdh{jj} = postinterp(fem,'VDragCoef',[dataPx{jj} dataPy{jj}].');
   dataAdh{jj} = dataAdh{jj}.';

   dataAdhFx{jj} = dataAdh{jj}.*dataUC1{jj};
   dataAdhFy{jj} = dataAdh{jj}.*dataUC2{jj};

   %Calculate the residue of the forword computed displacements.
   residueDisp{jj} = (dataUC1{jj}-dataU1{jj}).^2 + (dataUC2{jj}-dataU2{jj}).^2;

   %Calculate the residue of the regularized force.
   %residueBF{jj} = sigma*coef(:,jj).*coef(:,jj);
   residueBF{jj} = sigma*(coefMCFx(:,jj).*coefMCFx(:,jj) + ...
      coefMCFy(:,jj).*coefMCFy(:,jj)+coefAdh(:,jj).*coefAdh(:,jj));

   %The displacement on the edge to be compared with 'edgeU1' etc.
   for k = 1:numEdges
      [edgeUC1{jj,k} edgeUC2{jj,k}] = postinterp(fem,'u1','u2',edgeP{k});
   end

   %Save the identified body force calculated on the demonstration points.
   save([resultPath 'bfAdhId'],'bfDisplayPx','bfDisplayPy', ...
      'recMCF','recDispU','residueDisp','residueBF', ...
      'recADF','recAdhMap');

   %Save the computed displacement from the identified body force.
   save([resultPath 'dispAdhId'],'dataPx','dataPy','dataU1','dataU2', ...
      'dataUC1','dataUC2', 'dataMCFx','dataMCFy','dataAdh');
end

fprintf(1,'%f sec.\n',cputime-localStartTime);
