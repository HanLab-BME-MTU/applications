%Calculate the right-hand side vector of the system for identifying the body
% force. 
% This vector is given by substracting from the data the solution to the 
% elastic equation with zero body force and the inhomogeous Dirichlet boundary  
% condition where the boundary displacement is provided by the data.
% It can only run after running 'setupModel' and loading data.

fprintf(1,'Constructing the right-hand vector: \n'); 

startTime = cputime;

%iDispFieldDir = [reslDir filesep 'iDispField'];

%Set body force to be zero.
fp.BodyFx = {{'x' 'y'} {[] 0}};
fp.BodyFy = {{'x' 'y'} {[] 0}};

%Specify boundary condition.
%for k = 1:numEdges
%   fn.BndDispx{k} = 'bndDisp';
%   fn.BndDispy{k} = 'bndDisp';
%end
%fem = elModelUpdate(fem,'fn',fn,'fp',fp);

answer = input('Select time steps (0 for all):');
if isempty(answer) | answer == 0
   selTimeSteps = 1:numDTimePts;
else
   selTimeSteps = answer;
end

for ii = 1:length(selTimeSteps)
   jj = selTimeSteps(ii);

   localStartTime = cputime;
   fprintf(1,'   Time step %d ... ',jj);
   imgIndex = imgIndexOfDTimePts(jj);

   iDispFieldFileName = ['iDispField' sprintf(imgIndexForm,imgIndex) '.mat'];
   iDispFieldFile     = [iDispFieldDir filesep iDispFieldFileName];

   load(iDispFieldFile);
   numDP   = iDispField.numDP;
   iInMesh = iDispField.iInMesh;
   dataP   = iDispField.p(iInMesh,:).';
   edgD    = iDispField.edgD;

   if strcmp(dataToUse,'interp') == 1
      dataU1 = iDispField.v(iInMesh,1);
      dataU2 = iDispField.v(iInMesh,2);
   elseif strcmp(dataToUse,'smooth') == 1
      dataU1 = iDispField.sv(iInMesh,1);
      dataU2 = iDispField.sv(iInMesh,2);
   end

   if strcmp(isFieldBndFixed,'yes')
      femModelFile = [femModelDir filesep 'femModel' sprintf(imgIndexForm,0) '.mat'];
   else
      [femModelFile,femModelImgIndex] = getFemModelFile(femModelDir,imgIndex,imgIndexForm);
   end
   s = load(femModelFile);
   femModel = s.femModel;
   fem      = femModel.fem;
   numEdges = femModel.numEdges;

   for k = 1:numEdges
      fn.BndDispx{k} = 'bndDisp';
      fn.BndDispy{k} = 'bndDisp';
      fp.BndDispx{k} = {{'s'} {edgD(k).ppU1}};
      fp.BndDispy{k} = {{'s'} {edgD(k).ppU2}};
   end
   fem = elModelUpdate(fem,'fn',fn,'fp',fp);

   fem = elasticSolve(fem,[]);

   %'bspU1' and 'bspU2': Displacement with zero body force and boundary
   % displacement provided by the data.
   %[bspU1 bspU2] = postinterp(fem,'u1','u2',[dataPx{jj} dataPy{jj}].');
   [bspU1 bspU2] = postinterp(fem,'u1','u2',dataP);

   rightU = zeros(2*numDP,1);
   rightU(1:numDP) = dataU1-bspU1.';
   rightU(numDP+1:2*numDP) = dataU2-bspU2.';

   iDispField.bfRHU = rightU;

   %Save the modified 'iDispField' with calculated right hand U.
   save(iDispFieldFile,'iDispField');

   fprintf(1,'Done in %5.3f sec.\n',cputime-localStartTime);
   %For debugging, calculate the boundary displacements to compare with 
   % 'edge1U1' etc and calculate the body force to see if they are zeros.
   %for k = 1:numEdges
   %   [edgeUC1{jj,k} edgeUC2{jj,k}] = postinterp(fem,'u1','u2',edgeP{k});
   %end
   %[bspBFx bspBFy] = postinterp(fem,'f1','f2',[dataPx{jj} dataPy{jj}].');
end

fprintf(1,'Total time spent: %5.3f sec.\n', cputime-startTime);
