%Post assemble of the identified body force and some calculation
% for debugging purpose. It is run after 'solveLSBF'.

if strcmp(fwdOpComputed,'all') == 1
   fprintf(1,'Load the identified body force and displacements: \n');
   load([resultPath 'bfId']);
   load([resultPath 'dispId']);
   return;
end

fprintf(1,'Post assemble of the identified body force : ');

localStartTime = cputime;

load([modelPath 'femId']);
load([resultPath 'dispField']);
load([resultPath 'dataDisp']);
load([resultPath 'edgeDisp']);
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

%We assign scores to contraction and adhesion that identifies contraction site
% and adhesion sites with color.
scoreMCF = zeros(length(gridy),length(gridx));
scoreADF = scoreMCF;

%pp-form that interpolates 'scoreMCF' and 'scoreADF'.
ppMCF = cell(numTimeSteps,1);
ppADF = cell(numTimeSteps,1);

recBF       = zeros(length(bfDisplayPx),2,numTimeSteps);
residueDisp = cell(numTimeSteps,1);
residueBF   = cell(numTimeSteps,1);
recDispU    = zeros(length(bfDisplayPx),2,numTimeSteps);
dataUC1     = cell(numTimeSteps,1);
dataUC2     = cell(numTimeSteps,1);
dataBFx     = cell(numTimeSteps,1);
dataBFy     = cell(numTimeSteps,1);
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

   %Identify the location of adhesion and contration and assign intensity
   % scores that approximately separate the two. It is based on the angle
   % between the force and the displacement.
   %We first separate the two for forces on the display points.
   recBFLen = sqrt(recBFx.^2+recBFy.^2);
   recDispLen = sqrt(recDispU1.^2+recDispU2.^2);

   unitRecU1 = recDispU1./recDispLen;
   unitRecU2 = recDispU2./recDispLen;

   dotProdBFRecU = recBFx.*unitRecU1 + recBFy.*unitRecU2;
   mcfIndShow = find((dotProdBFRecU-recBFLen*cos(mcfAngle))>=0);
   adfIndShow = find((-dotProdBFRecU-recBFLen*cos(adfAngle))>=0);

   %We then do it on the grid points and assign intensity scores.
   %The force on grid points inside the identification region.
   gridXin = gridX(gridIn);
   gridYin = gridY(gridIn);
   [gridBFx,gridBFy] = postinterp(fem,'f1','f2', ...
      [gridX(gridIn) gridY(gridIn)].');

   %The displacements on grid points.
   [gridU1,gridU2] = postinterp(fem,'u1','u2', ...
      [gridX(gridIn) gridY(gridIn)].');

   %Compute the length of the body force.
   gridBFLen = sqrt(gridBFx.^2+gridBFy.^2);

   %Compute the length of the displacements and their average .
   gridDispLen = sqrt(gridU1.^2+gridU2.^2);
   smDispLen   = max(gridDispLen)*smDispThreshold;

   %Normalize the displacement.
   %unitGridU1 = zeros(size(gridU1));
   %unitGridU2 = zeros(size(gridU2));
   %unitGridU1(sigDispInd) = gridU1(sigDispInd)./gridDispLen(sigDispInd);
   %unitGridU2(sigDispInd) = gridU2(sigDispInd)./gridDispLen(sigDispInd);
   unitGridU1 = gridU1./gridDispLen;
   unitGridU2 = gridU2./gridDispLen;
   
   %Calculate the intensity scores for contraction and adhesion.
   % 'gridMCF' : Score for Myosin Contractile Force.
   % 'gridADF' : Score for Adhesion Dragging Force.
   %The sign of the dot product between the force and the displacement tells
   % if the angle between the force and the displacment is less or greater than
   % pi/2.
   %gridADF = -gridBFLen;
   %gridADF(sigDispInd) = gridBFx(sigDispInd).*unitGridU1(sigDispInd) + ...
   %   gridBFy(sigDispInd).*unitGridU2(sigDispInd);
   %First comput the dot product between the force and the displacements.
   dotProdBFU = gridBFx.*unitGridU1 + gridBFy.*unitGridU2;

   % mcfInd : Indices where we consider the total body force to be mainly 
   %          myosin contraction.
   % adfInd : Indices where we consider the total body force to be mainly 
   %          adhesion.
   mcfInd = find((dotProdBFU-gridBFLen*cos(mcfAngle))>=0);
   adfInd = find((-dotProdBFU-gridBFLen*cos(adfAngle))>=0);
   gridADF = zeros(size(gridBFLen));
   gridMCF = gridADF;

   %Significant myosin dragging force.
   gridMCF(mcfInd) = gridBFLen(mcfInd);
   
   %We only consider displacements that are above a threshold value. We call
   % these displacements the significant displacements.
   % 'smDispInd' : Index of the significant displacements.
   smDispInd = find(gridDispLen(adfInd)<smDispLen);
   %bigDispInd  = [1:length(gridDispLen)];
   %bigDispInd(smDispInd) = [];

   %Compute the significant dragging coefficient.
   gridADF(adfInd) = gridBFLen(adfInd)./gridDispLen(adfInd);
   %gridADF(smDispInd) = gridBFLen(smDispInd)./(avgDispLen*smDispThreshold);
   gridADF(adfInd(smDispInd)) = gridBFLen(adfInd(smDispInd))./smDispLen;

   %gridMCF = zeros(size(gridBFLen));
   %gridMCF(sigDispInd) = abs(gridBFx(sigDispInd).*unitGridU2(sigDispInd) - ...
   %   gridBFy(sigDispInd).*unitGridU1(sigDispInd));
   %gridMCF(mcfInd) = gridBFLen(mcfInd);

   %Interpolate 'gridMCF' and 'gridADF' to get the pp-form of the scores.
   scoreMCF(gridIn) = gridMCF./max(gridMCF);
   scoreADF(gridIn) = gridADF./max(gridADF);
   ppMCF{jj} = csape({gridy,gridx},scoreMCF);
   ppADF{jj} = csape({gridy,gridx},scoreADF);
   scoreMCF(:) = 0;
   scoreADF(:) = 0;

   %The displacements computed with the identified force and the identified
   % force at the data points. To be compared with 
   % 'dataU1' and 'dataU2'.
   [dataUC1{jj} dataUC2{jj}] = postinterp(fem,'u1','u2', ...
      [dataPx{jj} dataPy{jj}].');
   dataUC1{jj} = dataUC1{jj}.';
   dataUC2{jj} = dataUC2{jj}.';

   [dataBFx{jj} dataBFy{jj}] = postinterp(fem,'f1','f2', ...
      [dataPx{jj} dataPy{jj}].');
   dataBFx{jj} = dataBFx{jj}.';
   dataBFy{jj} = dataBFy{jj}.';

   %Calculate the residue of the forword computed displacements.
   residueDisp{jj} = (dataUC1{jj}-dataU1{jj}).^2 + (dataUC2{jj}-dataU2{jj}).^2;

   %Calculate the residue of the regularized force.
   %residueBF{jj} = sigma*coef(:,jj).*coef(:,jj);
   residueBF{jj} = sigma*(coefBFx(:,jj).*coefBFx(:,jj) + ...
      coefBFy(:,jj).*coefBFy(:,jj));

   %The displacement on the edge to be compared with 'edgeU1' etc.
   for k = 1:numEdges
      [edgeUC1{jj,k} edgeUC2{jj,k}] = postinterp(fem,'u1','u2',edgeP{k});
   end

   %Save the identified body force calculated on the demonstration points.
   save([resultPath 'bfId'],'bfDisplayPx','bfDisplayPy', ...
      'recBF','recDispU','residueDisp','residueBF', ...
      'gridXin','gridYin','gridBFx','gridBFy', ...
      'adfInd','mcfInd','ppMCF','ppADF');

   %Save the computed displacement from the identified body force.
   save([resultPath 'dispId'],'dataPx','dataPy','dataU1','dataU2', ...
      'dataUC1','dataUC2', 'dataBFx','dataBFy');
end

fprintf(1,'%f sec.\n',cputime-localStartTime);
