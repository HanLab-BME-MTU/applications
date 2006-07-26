%Calculate the right-hand side vector of the system for the boundary traction
% force. 
% This vector is given by substracting from the data the solution to the 
% elastic equation with given body force (already identified), 
% free boundary condition on the edge where the boundary 
% traction force is to be identified and the inhomogeous Dirichlet boundary 
% condition on the other edges where the boundary displacement is provided by 
% the data.
% It can only run after running 'setupModel' and loading data and after the
% body force is identified.

fprintf(1,'Constructing the right-hand vector:\n');

startTime = cputime;

%Set the body force.
%load([reslDir filesep 'bfId']);
%load([reslDir filesep 'coefBFId']);
fn.BodyFx = 'femBodyF';
fn.BodyFy = 'femBodyF';
fn.MyoDragFx = 0;
fn.MyoDragFy = 0;

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
   imgIndex = imgIndexOfDTimePts(jj);

   forceFieldFile = [forceFieldDir filesep 'forceField' ...
      sprintf(imgIndexForm,imgIndex) '.mat'];

   iDispFieldFileName = ['iDispField' sprintf(imgIndexForm,imgIndex) '.mat'];
   iDispFieldFile     = [iDispFieldDir filesep iDispFieldFileName];

   s = load(forceFieldFile);
   forceField = s.forceField;

   dataP  = forceField.p.';
   coefBF = forceField.coefBF;

   s = load(iDispFieldFile);
   iDispField = s.iDispField;

   numDP   = iDispField.numDP;
   iInMesh = iDispField.iInMesh;
   edgD    = iDispField.edgD;
   dataUC1 = iDispField.rv(:,1);
   dataUC2 = iDispField.rv(:,2);

   if strcmp(isFieldBndFixed,'yes')
      femModelFile = [femModelDir filesep 'femModel' sprintf(imgIndexForm,0) '.mat'];
   else
      [femModelFile,femModelImgIndex] = getFemModelFile(femModelDir,imgIndex,imgIndexForm);
   end
   s = load(femModelFile);
   femModel = s.femModel;
   fem      = femModel.fem;
   options  = femModel.options;
   fn       = femModel.fn;
   fp       = femModel.fp;
   fs       = femModel.fs;
   numEdges = femModel.numEdges;

   fp.BodyFx = {{'x' 'y'} {fs.fem coefBF(:,1)}};
   fp.BodyFy = {{'x' 'y'} {fs.fem coefBF(:,2)}};
   for k = 1:numEdges
      BCTypes{k}     = 'Dirichlet';
      fn.BndDispx{k} = 'bndDisp';
      fn.BndDispy{k} = 'bndDisp';
      fp.BndDispx{k} = {{'s'} {edgD(k).ppU1}};
      fp.BndDispy{k} = {{'s'} {edgD(k).ppU2}};

      fn.BndTracFx{k} = 0;
      fn.BndTracFy{k} = 0;
   end

   procStr = '';
   rightU  = cell(numEdges,1);
   for k = 1:numEdges
      BCTypes{k} = 'Neumann';
      options    = elOptionsSet(options,'BCType',BCTypes);

      fem = elModelUpdate(fem,'options',options,'fn',fn,'fp',fp);
      fem = elasticSolve(fem,[]);

      %Substract the solution from the data to get the right hand vector.
      [bspU1 bspU2] = postinterp(fem,'u1','u2',dataP);

      rightU{k} = zeros(2*numDP,1);
      rightU{k}(1:numDP) = dataUC1-bspU1.';
      rightU{k}(numDP+1:2*numDP) = dataUC2-bspU2.';

      BCTypes{k} = 'Dirichlet';

      %For debugging, calculate the body force. It is supposed to equal 'recBFx' and 'recBFy'.
      [tfRecBFx,tfRecBFy] = postinterp(fem,'f1','f2',dataP);

      %For debugging, calculate the boundary displacements to compare with 'edgeU1'
      % and 'edgeU2'.
      %for j = 1:numEdges
      %   [edgeUC1{j} edgeUC2{j}] = postinterp(fem,'u1','u2',edgeP{j});
      %end

      for ik = 1:length(procStr)
         fprintf(1,'\b');
      end
      procStr = sprintf('   Edge %d (out of %d) finished in %5.3f sec.', ...
         k,numEdges,cputime-localStartTime);
      fprintf(1,procStr);
   end
   fprintf(1,'\n');

   iDispField.tfRHU = rightU;

   %Save the modified 'iDispField' with calculated right hand U.
   save(iDispFieldFile,'iDispField');
end

fprintf(1,'Total time spent : %5.3f sec.\n', cputime-startTime);
