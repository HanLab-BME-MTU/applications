%This script file computes the matrix approximation to the forward linear
% operator for the body force. It can only run after running 'setupModel' and 
% loading data.
% For the function space of the boundary traction force, we use B-spline since
% it it one dimention.

%edgBrksTF : The breaks on each edge for the construction of B-spline.
%edgKnotsT : FThe knot sequence on each edge built from the breaks.
%dimTF     : Dimention of boundary Traction Force.
%bspTF     : The sp-form of the boundary Traction Force.
%coefTF    : The coefficient for the construction of the boundary Traction 
%            Force.
%A         : The matrix respresentation of the forward operator for the
%            boundary Traction Force.

fprintf(1,['\nCalculating the forward operator for the boundary ' ...
   'traction force :\n']);

startTime = cputime;

if strcmp(tfFwdOpComputed,'none') == 1
   localStartTime = cputime;

   %Load 'femModel' that contains the field geometry, fem model and domain force function space.
   [modelFileList,modelImgIndex] = getFemModelFile(femModelDir);
   if isempty(modelFileList)
      fprintf(1,'Fem model has not been set up and saved yet. Run setupModel first.\n');
      return;
   end

   answer = input('Select Model fields (0 for all):');
   if isempty(answer) | answer == 0
      selFields = 1:length(modelFileList);
   else
      selFields = answer;
   end

   for ii = 1:length(selFields)
      jj = selFields(ii);

      localStartTime = cputime;

      imgIndex = modelImgIndex(jj);

      procStr = '';
      fprintf(1,'  Field %d: ',jj);

      s = load(modelFileList{jj});
      femModel = s.femModel;
      numEdges = femModel.numEdges;
      edge     = femModel.edge;
      fem      = femModel.fem;
      fn       = femModel.fn;
      fp       = femModel.fp;
      options  = femModel.options;

      DTDir = ['DT' sprintf(imgIndexForm,imgIndex)]; 
      if ~isdir([femSolBasisTFDir filesep DTDir])
         [success,msg,msgId] = mkdir(femSolBasisTFDir,DTDir);
         if ~success
            error('Trouble making directory for ''femSolBasisTF''.');
         end
      end
      solTFIdFile = [femSolBasisTFDir filesep DTDir filesep 'solId'];
      %fsBndFile   = [femSolBasisTFDir filesep DTDir filesep 'fsBnd.mat'];

      for k = 1:numEdges
         numEdgBrksTF  = floor(edge(k).arcLen/edgBrkDistTF);

         fsBnd(k).edgBrks  = linspace(0,edge(k).arcLen,numEdgBrksTF);
         fsBnd(k).edgKnots = augknt(fsBnd(k).edgBrks,bspOrderTF);
         fsBnd(k).dim      = length(fsBnd(k).edgKnots)-bspOrderTF;

         %coefTF{k} = zeros(1,fsBnd(k).dim);
         %bspTF{k}  = cell(dimTF(k),1);
      end
      femModel.fsBnd = fsBnd;
      save(modelFileList{jj},'femModel');

      BCTypes = options.BCType;
      for k = 1:numEdges
         %Construct the sequence of breaks according to 'EdgBrkDistTF', the
         % distance between the breaks.
         sol{k} = cell(2*fsBnd(k).dim,1);

         BCTypes{k} = 'Neumann';
         options    = elOptionsSet(options,'BCType',BCTypes);

         fem = elModelUpdate(fem,'options',options);

         coefTF = zeros(1,fsBnd(k).dim);
         for j = 1:fsBnd(k).dim
            coefTF(j) = 1;
            bspTF = spmak(fsBnd(k).edgKnots,coefTF);
            coefTF(j) = 0;

            fp.BndTracFx{k} = {{'s'} {bspTF}};
            fp.BndTracFy{k} = {{'s'} {0}};
            fem = elModelUpdate(fem,'fp',fp);
            fem = elasticSolve(fem,[]);
            sol{k}{2*j-1} = fem.sol;

            %'bspU1' and 'bspU2' is used to temporarily store the solution at the
            % data points.
            %[bspU1 bspU2] = postinterp(fem,'u1','u2',[dataPx{jj} dataPy{jj}].');
            %A{k}(:,1,j,1)  = bspU1.';
            %A{k}(:,2,j,1)  = bspU2.';

            fp.BndTracFx{k} = {{'s'} {0}};
            fp.BndTracFy{k} = {{'s'} {bspTF}};
            fem = elModelUpdate(fem,'fp',fp);
            fem = elasticSolve(fem,[]);
            sol{k}{2*j} = fem.sol;

            %'bspU1' and 'bspU2' is used to temporarily store the solution at the
            % data points.
            %[bspU1 bspU2] = postinterp(fem,'u1','u2',[dataPx{jj} dataPy{jj}].');
            %A{k}(:,1,j,2)  = bspU1.';
            %A{k}(:,2,j,2)  = bspU2.';
         end
         for kk = 1:length(procStr)
            fprintf(1,'\b');
         end
         procStr = sprintf('   Edge %d (out of %d) finished in %5.3f sec.', ...
            k,numEdges,cputime-localStartTime);
         fprintf(1,procStr);

         BCTypes{k} = 'Dirichlet';
      end
      fprintf(1,'\n');

      %save([reslDir filesep 'AtfId'],'A');
      %save([mechDir filesep 'solTFId'],'sol');
      save(solTFIdFile,'sol');
   end
end

if ~strcmp(tfFwdOpComputed,'A') || ~strcmp(tfFwdOpComputed,'all')
   answer = input('Select time steps (0 for all):');
   if isempty(answer) | answer == 0
      selTimeSteps = 1:numDTimePts;
   else
      selTimeSteps = answer;
   end

   fprintf(1,'Constructing the matrix A for boundary force:\n');
   clear A;
   for ii = 1:length(selTimeSteps)
      jj = selTimeSteps(ii);

      procStr = '';
      fprintf(1,'  Time step %d: ', jj);
      localStartTime = cputime;

      imgIndex = imgIndexOfDTimePts(jj);
      iDispFieldFileName = ['iDispField' sprintf(imgIndexForm,imgIndex) '.mat'];
      iDispFieldFile     = [iDispFieldDir filesep iDispFieldFileName];
      s = load(iDispFieldFile);
      iDispField = s.iDispField;

      numDP = iDispField.numDP;

      if strcmp(isFieldBndFixed,'yes')
         femModelFile = [femModelDir filesep 'femModel' sprintf(imgIndexForm,0) '.mat'];
         DTDir = ['DT' sprintf(imgIndexForm,0)];
      else
         [femModelFile,femModelImgIndex] = getFemModelFile(femModelDir,imgIndex,imgIndexForm);
         DTDir = ['DT' sprintf(imgIndexForm,femModelImgIndex)]; 
      end
      s = load(femModelFile);
      femModel = s.femModel;

      fsBnd = femModel.fsBnd;
      fem   = femModel.fem;

      solTFIdFile = [femSolBasisTFDir filesep DTDir filesep 'solId'];
      %fsBndFile   = [femSolBasisTFDir filesep DTDir filesep 'fsBnd.mat'];

      s = load(solTFIdFile,'sol');
      sol = s.sol;
      for k = 1:numEdges
         %Construct the sequence of breaks according to 'EdgBrkDistTF', the
         % distance between the breaks.
         A{k}    = zeros(numDP,2,fsBnd(k).dim,2);

         dataP = iDispField.p(iDispField.iInMesh,:).';
         for j = 1:fsBnd(k).dim
            %'bspU1' and 'bspU2' is used to temporarily store the solution at the
            % data points.
            fem.sol = sol{k}{2*j-1};
            [bspU1 bspU2] = postinterp(fem,'u1','u2',dataP);
            A{k}(:,1,j,1) = bspU1.';
            A{k}(:,2,j,1) = bspU2.';

            fem.sol = sol{k}{2*j};
            %[bspU1 bspU2] = postinterp(fem,'u1','u2',[dataPx{jj} dataPy{jj}].');
            [bspU1 bspU2] = postinterp(fem,'u1','u2',dataP);
            A{k}(:,1,j,2)  = bspU1.';
            A{k}(:,2,j,2)  = bspU2.';
         end
         for kk = 1:length(procStr)
            fprintf(1,'\b');
         end
         procStr = sprintf('   Edge %d (out of %d) finished in %5.3f sec.', ...
            k,numEdges,cputime-localStartTime);
         fprintf(1,procStr);
      end
      fprintf(1,'\n');

      fwdMapTFFile = [fwdMapTFDir filesep 'A' sprintf(imgIndexForm,imgIndex) '.mat'];
      save(fwdMapTFFile,'A');
   end
end


