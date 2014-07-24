%Post assemble of the identified boundary traction force and some calculation
% for debugging purpose. It is run after 'solveLSTF'.

fprintf(1,'Post assemble of the identified boundary traciton force:\n');

startTime = cputime;

ans = input('Select time steps (0 for all):');
if isempty(ans) | ans == 0
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
   indexStr = sprintf(imgIndexForm,imgIndex);

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
   options  = femModel.options;

   forceFieldFile = [forceFieldDir filesep 'forceField' indexStr '.mat'];
   s = load(forceFieldFile);
   forceField = s.forceField;

   iDispFieldFileName = ['iDispField' indexStr '.mat'];
   iDispFieldFile     = [iDispFieldDir filesep iDispFieldFileName];
   s = load(iDispFieldFile);
   iDispField = s.iDispField;

   iDispField.edgDispIntvLen = edgDispIntvLen;

   edgD = iDispField.edgD;

   bndfVecFile  = [reslDir filesep 'bndfVecF' filesep 'bndfVecF' indexStr '.mat'];
   bndfMagFile  = [reslDir filesep 'bndfMag' filesep 'bndfMag' indexStr '.mat'];
   
   if ~exist([reslDir filesep 'bndfVecF'],'dir')
      [success msg msgID] = mkdir(reslDir,'bndfVecF');
      if ~success
         error('Trouble making directory ''bndfVecF''.');
      end
   end

   if ~exist([reslDir filesep 'bndfMag'],'dir')
      [success msg msgID] = mkdir(reslDir,'bndfMag');
      if ~success
         error('Trouble making directory ''bndfMag''.');
      end
   end
   
   bndfVec = [];
   bndfMag = [];
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
      %edgCornerSegCut = 4;
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

      bndfVec = [bndfVec; [bndF.p(2:-1:1,:); bndF.p(2:-1:1,:)+[bndF.fy; bndF.fx]].'];
      bndfMag = [bndfMag; [imgIndex*ones(1,size(bndF.p,2)); bndF.p(2:-1:1,:); sqrt(bndF.fx.^2+bndF.fy.^2)].'];
      forceField.bndF = [forceField.bndF bndF];
   end
   save(bndfVecFile,'bndfVec');
   save(bndfMagFile,'bndfMag');
   save(forceFieldFile,'forceField');

   %Calculating the flow field given the reconstructed domain force and reconstructed boundary force
   % on one edge at a time.
   fem = femModel.fem;
   fs  = femModel.fs;
   fn  = femModel.fn;
   fp  = femModel.fp;

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

   BCTypes = options.BCType;
   for k = 1:numEdges
      BCTypes{k} = 'Neumann';
      options    = elOptionsSet(options,'BCType',BCTypes);

      fem = elModelUpdate(fem,'options',options,'fn',fn,'fp',fp);
      fem = elasticSolve(fem,[]);

      %The displacements computed with the identified force. To be compared with 
      % 'dataU1' and 'dataU2'.

      [UC1 UC2] = postinterp(fem,'u1','u2',forceField.p.');
      iDispField.bndF_rv{k} = [UC1;UC2].';

      [edgeUC1 edgeUC2] = postinterp(fem,'u1','u2',edge(k).bndP);
      iDispField.rEdgD{k} = [edgeUC1;edgeUC2].';

      %The boundary traction force on the edge to be compared with 'TFxR' etc.
      %[edgeFx edgeFy] = postinterp(fem,'g1','g2',edge(k).bndP);

      BCTypes{k} = 'Dirichlet';
   end
   save(iDispFieldFile,'iDispField');

   fprintf(1,'Done in %5.3f sec.\n',cputime-localStartTime);
end

fprintf(1,'Total time spent: %5.3f sec.\n',cputime-startTime);
