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

%Get the grid points that are inside the identification region. They will be
% used to calculate the intensity score for contraction and adhesion.
[is,pe] = postinterp(fem,[gridX gridY].');
gridIn = 1:length(gridX);
gridIn(pe) = [];

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
   %The force on grid points inside the identification region.
   [gridBFx,gridBFy] = postinterp(fem,'f1','f2', ...
      [gridX(gridIn) gridY(gridIn)].');

   %The displacements on grid points.
   [gridU1,gridU2] = postinterp(fem,'u1','u2', ...
      [gridX(gridIn) gridY(gridIn)].');

   %Compute the length of the body force.
   gridBFLen = sqrt(gridBFx.^2+gridBFy.^2);

   %Compute the length of the displacements and their average .
   gridDispLen = sqrt(gridU1.^2+gridU2.^2);
   avgDispLen = sum(gridDispLen(:))/length(gridDispLen(:));

   %We only consider displacements that are above a threshold value. We call
   % these displacements the significant displacements.
   % 'sigDispInd' : Index of the significant displacements.
   sigDispInd = find(gridDispLen>=avgDispLen*sigDispThreshold);
   smDispInd  = [1:length(gridDispLen)];
   smDispInd(sigDispInd) = [];

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
   gridADF = gridBFx.*unitGridU1 + gridBFy.*unitGridU2;

   % mcfInd : Indices where we assign the total body force to myosin
   % contraction and zero to adhesion.
   mcfInd = find(gridADF>=0);
   gridADF(mcfInd) = 0;
   gridADF = abs(gridADF);
   
   %Compute the dragging coefficient.
   gridADF(smDispInd)  = gridADF(smDispInd)./(avgDispLen*sigDispThreshold);
   gridADF(sigDispInd) = gridADF(sigDispInd)./gridDispLen(sigDispInd);

   gridMCF = zeros(size(gridBFLen));
   gridMCF(sigDispInd) = abs(gridBFx(sigDispInd).*unitGridU2(sigDispInd) - ...
      gridBFy(sigDispInd).*unitGridU1(sigDispInd));
   gridMCF(mcfInd) = gridBFLen(mcfInd);

   %Interpolate 'gridMCF' and 'gridADF' to get the pp-form of the scores.
   scoreMCF(gridIn) = gridMCF;
   scoreADF(gridIn) = gridADF;
   ppMCF{jj} = csape({gridy,gridx},scoreMCF);
   ppADF{jj} = csape({gridy,gridx},scoreADF);
   scoreMCF(:) = 0;
   scoreADF(:) = 0;

   %The displacements computed with the identified force. To be compared with 
   % 'dataU1' and 'dataU2'.
   [dataUC1{jj} dataUC2{jj}] = postinterp(fem,'u1','u2', ...
      [dataPx{jj} dataPy{jj}].');
   dataUC1{jj} = dataUC1{jj}.';
   dataUC2{jj} = dataUC2{jj}.';

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
      'recBF','recDispU','residueDisp','residueBF','ppMCF','ppADF');

   %Save the computed displacement from the identified body force.
   save([resultPath 'dispId'],'dataPx','dataPy','dataUC1','dataUC2');
end

fprintf(1,'%f sec.\n',cputime-localStartTime);
