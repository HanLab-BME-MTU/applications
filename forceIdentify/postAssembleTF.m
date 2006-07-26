%Post assemble of the identified boundary traction force and some calculation
% for debugging purpose. It is run after 'solveLSTF'.

fprintf(1,'Post assemble of the identified boundary traciton force:\n');

startTime = cputime;

ans = input('Select time steps (0 for all):');
if isempty(ans) || ans == 0
   selTimeSteps = 1:numDTimePts;
else
   selTimeSteps = ans;
end

for ii = 1:length(selTimeSteps)
   jj = selTimeSteps(ii);

   localStartTime = cputime;

   fprintf(1,'   Time step %d ... ',jj);

   %Get one cell image
   imgIndex = imgIndexOfDTimePts(jj);

   if strcmp(isFieldBndFixed,'yes')
      femModelFile = [femModelDir filesep 'femModel' sprintf(imgIndexForm,0) '.mat'];
   else
      [femModelFile,femModelImgIndex] = getFemModelFile(femModelDir,imgIndex,imgIndexForm);
   end
   s = load(femModelFile);
   femModel = s.femModel;
   edge     = femModel.edge;
   numEdges = femModel.numEdges;
   fem      = femModel.fem;
   fn       = femModel.fn;
   fp       = femModel.fp;
   fsBnd    = femModel.fsBnd;

   forceFieldFile = [forceFieldDir filesep 'forceField' ...
      sprintf(imgIndexForm,imgIndex) '.mat'];
   s = load(forceFieldFile);
   forceField = s.forceField;

   iDispFieldFileName = ['iDispField' sprintf(imgIndexForm,imgIndex) '.mat'];
   iDispFieldFile     = [iDispFieldDir filesep iDispFieldFileName];
   s = load(iDispFieldFile);
   iDispField = s.iDispField;

   iDispField.edgDispIntvLen = edgDispIntvLen;

   edgD = iDispField.edgD;

   forceField.bndF = [];
   for k = 1:numEdges
      clear bndF;

      numEdgDispPt = floor(edge(k).arcLen/edgDispIntvLen);
      bndF.numDispPt = numEdgDispPt;
      bndF.s = linspace(0,edge(k).arcLen,numEdgDispPt);
      bndF.p = [fnval(edge(k).ppX,bndF.s); ...
                fnval(edge(k).ppY,bndF.s)];

      coefTFx = forceField.coefTF{k}(1:fsBnd(k).dim).';
      coefTFy = forceField.coefTF{k}(fsBnd(k).dim+1:end).';

      %The corner between two edges poses a singular point for boundary force reconstruction.
      % We cut off e.g. 2 segments.
      edgCornerSegCut = 4;
      coefTFx(1:edgCornerSegCut) = 0;
      coefTFy(1:edgCornerSegCut) = 0;
      coefTFx(end-edgCornerSegCut:end) = 0;
      coefTFy(end-edgCornerSegCut:end) = 0;
      bndF.spx = spmak(fsBnd(k).edgKnots,coefTFx);
      bndF.spy = spmak(fsBnd(k).edgKnots,coefTFy);

      bndF.fx = fnval(bndF.spx,bndF.s);
      bndF.fy = fnval(bndF.spy,bndF.s);

      bndF.amp = sqrt(bndF.fx.^2+bndF.fy.^2);
      bndF.phs = angle(bndF.fx+i*bndF.fy);

      forceField.bndF = [forceField.bndF bndF];
   end
   save(forceFieldFile,'forceField');

   if strcmp(debugMode,'on')
      %For debugging : Solve the elastic equation with the identified force and
      % compare the computed displacement with the measured displacement data.
      fem = femModel.fem;
      fs  = femModel.fs;

      %Set the body force.
      coefBFx = forceField.coefBF(:,1);
      coefBFy = forceField.coefBF(:,2);

      fn.BodyFx = 'femBodyF';
      fn.BodyFy = 'femBodyF';
      fp.BodyFx = {{'x' 'y'} {fs.fem coefBFx}};
      fp.BodyFy = {{'x' 'y'} {fs.fem coefBFy}};

      %Specify boundary condition.
      for k = 1:numEdges
         fn.BndDispx{k} = 'bndDisp';
         fn.BndDispy{k} = 'bndDisp';
         fp.BndDispx{k} = {{'s'} {edgD(k).ppU1}};
         fp.BndDispy{k} = {{'s'} {edgD(k).ppU2}};

         bndF = forceField.bndF(k);
         fn.BndTracFx{k} = 'spBndTF';
         fn.BndTracFy{k} = 'spBndTF';
         fp.BndTracFx{k} = {{'s'} {bndF.spx}};
         fp.BndTracFy{k} = {{'s'} {bndF.spy}};
      end

      for k = 1:numEdges
         BCTypes{k} = 'Neumann';
         options    = elOptionsSet(options,'BCType',BCTypes);

         fem = elModelUpdate(fem,'options',options,'fn',fn,'fp',fp);
         fem = elasticSolve(fem,[]);

         %The displacements computed with the identified force. To be compared with 
         % 'dataU1' and 'dataU2'.

         [edgeUC1 edgeUC2] = postinterp(fem,'u1','u2',edge(k).bndP);
         iDispField.rEdgD{k} = [edgeUC1;edgeUC2].';

         %The boundary traction force on the edge to be compared with 'TFxR' etc.
         %[edgeFx edgeFy] = postinterp(fem,'g1','g2',edge(k).bndP);

         BCTypes{k} = 'Dirichlet';
      end
      save(iDispFieldFile,'iDispField');
   end

   fprintf(1,'Done in %5.3f sec.\n',cputime-localStartTime);
end

fprintf(1,'Total time spent: %5.3f sec.\n',cputime-startTime);
