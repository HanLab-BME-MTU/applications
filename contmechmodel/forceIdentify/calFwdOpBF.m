%This script file computes the matrix approximation to the forward linear
% operator for the body force. It can only run after running 'setupModel' and 
% 'calFwdOpBF' or loading data.

fprintf(1,'Calculating the forward operator for the Body Force :\n');

startTime = cputime;

if strcmp(bfFwdOpComputed,'none') == 1
   %Load 'femModel' that contains the field geometry, fem model and domain force function space.
   [modelFileList,modelImgIndex] = getFemModelFile(femModelDir);
   if isempty(modelFileList)
      fprintf(1,'Fem model has not been set up and saved yet. Run setupModel first.\n');
      return;
   end

   answer = input('Select Model fields (0 for all):');
   if isempty(answer)
      selFields = 1:length(modelFileList);
   elseif length(answer) == 1 && answer == 0
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
      fem      = femModel.fem;
      fs       = femModel.fs;

      %Step 2: Construct the matrix approximation to the forward operator.
      %We use each basis function as the body force to solve our
      % continuum mechanics system. The solution gives us each column of the matrix.
      dimFS     = fs.dimFS;
      dimBF     = fs.dimBF;
      indDomDOF = fs.indDomDOF;
      coefFS    = zeros(dimFS,1);
      solFileIndexForm = sprintf('%%.%dd',length(num2str(2*dimBF)));
      if rem(2*dimBF,numBSolsPerFile) == 0
         numSolFiles = 2*dimBF/numBSolsPerFile;
      else
         numSolFiles = ceil(2*dimBF/numBSolsPerFile);
      end

      %fprintf(1,'\n  Solving the elastic equation for each basis function :\n');
      %To store the solutions corresponding to each basis.
      %sol = cell(2*dimBF,1); 
      sol = cell(numBSolsPerFile,1); 

      DTDir = ['DT' sprintf(imgIndexForm,imgIndex)]; 

      %col : goes from 1 to 2*dimBF
      %ll  : goes from 1 to numBSolsPerFile.
      col = 0;
      ll  = 0;
      startCol = 1;
      j = 1;
      while j <= dimBF
         %Check if the basis solution has been calculated and saved before.
         if isdir([femSolBasisBFDir filesep DTDir]) && strcmp(overwriteBSolBF,'no')
            endCol = min(2*dimBF, startCol+numBSolsPerFile-1);
            solBFIdFile = [femSolBasisBFDir filesep DTDir filesep 'solId' ...
               sprintf(solFileIndexForm,startCol) '_' ...
               sprintf(solFileIndexForm,endCol) '.mat'];

            if exist(solBFIdFile,'file') == 2
               j        = j+(endCol-startCol+1)/2;
               startCol = endCol+1;
               ll       = 0;
               col      = endCol;

               continue;
            end
         end

         %'k' is the index of 'coefFS' whose corresponding basis function is zero
         % on the boundary.
         k = indDomDOF(j);
         coefFS(k) = 1;
         fp.BodyFx = {{'x' 'y'} {fs.fem coefFS}};
         fp.BodyFy = {{'x' 'y'} {[] 0}};
         fem = elModelUpdate(fem,'fp',fp);
         fem = elasticSolve(fem,[]);
         sol{ll+1} = fem.sol;

         fp.BodyFx = {{'x' 'y'} {[] 0}};
         fp.BodyFy = {{'x' 'y'} {fs.fem coefFS}};
         fem = elModelUpdate(fem,'fp',fp);
         fem = elasticSolve(fem,[]);
         sol{ll+2} = fem.sol;

         coefFS(k) = 0;

         col = col+2;
         ll  = ll+2;
         if rem(col,10) == 0
            for kk = 1:length(procStr)
               fprintf(1,'\b');
            end
            procStr = sprintf('%d basis functions out of %d finished in %5.3f sec ...', ...
               col,2*dimBF,cputime-localStartTime);
            fprintf(1,procStr);
         end

         if ll == numBSolsPerFile || col == 2*dimBF
            endCol = col;
            if ~isdir([femSolBasisBFDir filesep DTDir])
               success = mkdir(femSolBasisBFDir,DTDir);
               if ~success
                  error('Trouble making DT directory.');
               end
            end
            solBFIdFile = [femSolBasisBFDir filesep DTDir filesep 'solId' ...
               sprintf(solFileIndexForm,startCol) '_' ...
               sprintf(solFileIndexForm,endCol) '.mat'];
            save(solBFIdFile,'sol');

            %Restart the counting of 'll'.
            ll  = 0;
            sol = cell(numBSolsPerFile,1);
            startCol = endCol+1;
         end
         j = j+1;
      end %Of j<= dimBF.

      fprintf(1,'\n');
      %elseif strcmp(bfFwdOpComputed,'fem') == 1
      %   fprintf(1,['  Loading the solution to the elastic equation ' ...
      %      'for each basis function.\n']);
      %   load([mechDir filesep  'solBFId'],'sol');
   end
end

%Evaluate 'femSolBasisBF' at data points to get columns of the matrix representation of the
% forward map.
%if strcmp(bfFwdOpComputed,'A') == 1 | strcmp(bfFwdOpComputed,'all')
%   fprintf(1,'  Loading the matrix A for the Body Force.\n');
%   load([reslDir filesep  'AbfId'],'A');
%else
if ~strcmp(bfFwdOpComputed,'A') && ~strcmp(bfFwdOpComputed,'all')
   if strcmp(isFieldBndFixed,'yes') && strcmp(isDataPosFixed,'yes')
      %We only need to calculate 'A' once for all time points in this case.
      selTimeSteps = 0;
   else
      answer = input('Select time steps (0 for all):');
      if isempty(answer) | answer == 0
         selTimeSteps = 1:numDTimePts;
      else
         selTimeSteps = answer;
      end
   end

   fprintf(1,'Constructing the matrix A for domain force:\n');
   clear A;
   for ii = 1:length(selTimeSteps)
      jj = selTimeSteps(ii);

      localStartTime = cputime;
      procStr = '';
      if jj == 0
         fprintf(1,'  Both domain and data positions are fixed: for all time steps.\n');
         imgIndex = imgIndexOfDTimePts(1);
      else
         fprintf(1,'  Time step %d: ', jj);
         imgIndex = imgIndexOfDTimePts(jj);
      end

      %Load the interpolated displacement field.
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
      fem      = femModel.fem;
      fs       = femModel.fs;

      %Step 2: Construct the matrix approximation to the forward operator.
      %We use each basis function as the body force to solve our
      % continuum mechanics system. The solution gives us each column of the matrix.
      dimFS     = fs.dimFS;
      dimBF     = fs.dimBF;
      indDomDOF = fs.indDomDOF;
      coefFS    = zeros(dimFS,1);
      solFileIndexForm = sprintf('%%.%dd',length(num2str(2*dimBF)));
      if rem(2*dimBF,numBSolsPerFile) == 0
         numSolFiles = 2*dimBF/numBSolsPerFile;
      else
         numSolFiles = ceil(2*dimBF/numBSolsPerFile);
      end

      A = zeros(numDP,2,dimBF,2);

      startCol = 1;
      ll  = 0;
      for j = 1:numSolFiles
         endCol = min(2*dimBF,j*numBSolsPerFile);

         solBFIdFile = [femSolBasisBFDir filesep DTDir filesep 'solId' ...
            sprintf(solFileIndexForm,startCol) '_' ...
            sprintf(solFileIndexForm,endCol) '.mat'];

         s = load(solBFIdFile,'sol');
         sol = s.sol;
         %solBFIdFile = ['solBF' filesep 'solId' ...
         %   sprintf(solFileIndexForm, startCol) '_' ...
         %   sprintf(solFileIndexForm,endCol)];
         %load([mechDir filesep  solBFIdFile], 'sol');
         for col = 1:2:endCol-startCol
            ll = ll+1;

            %'bspU1' and 'bspU2' is used to temporarily store the solution at 
            % the data points.
            fem.sol = sol{col};
            %[bspU1 bspU2] = postinterp(fem,'u1','u2',[dataPx{jj} dataPy{jj}].');
            dataP = iDispField.p(iDispField.iInMesh,:).';

            [bspU1 bspU2] = postinterp(fem,'u1','u2',dataP);
            A(:,1,ll,1)  = bspU1.';
            A(:,2,ll,1)  = bspU2.';

            fem.sol = sol{col+1};
            %[bspU1 bspU2] = postinterp(fem,'u1','u2',[dataPx{jj} dataPy{jj}].');
            [bspU1 bspU2] = postinterp(fem,'u1','u2',dataP);
            A(:,1,ll,2)  = bspU1.';
            A(:,2,ll,2)  = bspU2.';

            if rem(2*ll,100) == 0
               for kk = 1:length(procStr)
                  fprintf(1,'\b');
               end
               procStr = sprintf('columns %d out of %d finished in %5.3f sec ...', ...
                  2*ll,2*dimBF,cputime-localStartTime);
               fprintf(1,procStr);
            end
         end
         startCol = endCol+1;
      end
      fprintf(1,'\n');

      %Reshape 'A';
      A = reshape(A,2*numDP,2*dimBF);

      if jj == 0
         fwdMapBFFile = [fwdMapBFDir filesep 'A' sprintf(imgIndexForm,0) '.mat'];
      else
         fwdMapBFFile = [fwdMapBFDir filesep 'A' sprintf(imgIndexForm,imgIndex) '.mat'];
      end
      save(fwdMapBFFile,'A');
      %save([reslDir filesep  'AbfId'],'A');
   end
end

fprintf(1,'Total time spent: %5.3f sec.\n', cputime-startTime);

%Some old code.
%   if exist('fsMesh')==1 && strcmp(fsMesh,'femMesh')
%      curve = geomspline([bfDomPGx bfDomPGy].', ...
%      'splinemethod','centripetal','closed','auto');
%      geom  = geomcoerce('solid',{curve});
%
%      if exist('bfConstraint') == 1 &% strcmp(bfConstraint,'adhLocation') == 1
%         if strcmp(adhGeomType,'node') == 1
%            fs.geom = geomcoerce('solid',{curve adhNodes{:}});
%         else
%            fs.geom = geom+adhGeom;
%         end
%
%         %Get the number of subdomains.
%         numSubDoms = geominfo(fs.geom,'Out',{'nmr'});
%
%         %Specify a smaller enought meshing size for adhesion sites. Note:
%         % Subdomain No. 1 should be the main domain surrounding the adhesion
%         % sites.
%         subHmax = zeros(2,numSubDoms-1);
%         for k = 2:numSubDoms
%            subHmax(:,k-1) = [k;4*sqrt(smAdhThreshold)];
%         end
%
%         %Get the global 'hmax' if it is specified in 'fsMeshPar'.
%         %The position of the 'hmax' property.
%         if exist('fsMeshPar')
%            hmaxId = find(strcmp(fsMeshPar,'hmax'));
%         else
%            fsMeshPar = {};
%            hmaxId    = [];
%         end
%
%         if isempty(hmaxId)
%            %A gloabl 'hmax' is not specified in 'fsMeshPar'.
%            fsMeshPar{end+1} = 'hmax';
%            fsMeshPar{end+1} = {3*sqrt(smAdhThreshold) [] [] subHmax};
%            %fsMeshPar{end+1} = 2*sqrt(smAdhThreshold);
%         else
%            %Merge the global 'hmax' with the adhesion 'hmax'.
%            fsMeshPar{hmaxId+1} = {fsMeshPar{hmaxId+1} [] [] subHmax};
%         end
%      else
%         fs.geom = geom;
%      end
%
%      if exist('fsMeshPar')
%         fs.mesh = meshinit(fs,fsMeshPar{:});
%      else
%         fs.mesh = meshinit(fs);
%      end
%   else
%      curvL = [bfDomPGx(bfDomPGVI(1):bfDomPGVI(2)) ...
%               bfDomPGy(bfDomPGVI(1):bfDomPGVI(2))].';
%      curvT = [bfDomPGx(bfDomPGVI(2):bfDomPGVI(3)) ...
%               bfDomPGy(bfDomPGVI(2):bfDomPGVI(3))].';
%      curvR = [bfDomPGx(bfDomPGVI(3):bfDomPGVI(4)) ...
%               bfDomPGy(bfDomPGVI(3):bfDomPGVI(4))].';
%      curvB = [bfDomPGx(bfDomPGVI(4):end) bfDomPGy(bfDomPGVI(4):end)].';
%      fs.mesh  = rectmesh(curvL,curvT,curvR,curvB,bfDomHin,bfDomVin);
%      fs.geom  = fs.mesh;
%   end

   %For backward compatibility, we load the one file of solutions and break it
   % into multiple files
   %if strcmp(bfFwdOpComputed,'fem') == 1 && ~isdir([mechDir filesep 'solBF']) ...
   %   && exist([mechDir filesep 'solBFId.mat'],'file') == 2
   %   S = load([mechDir filesep  'solBFId'],'sol');
   %   oldSol = S.sol;
   %
   %   [success,msg,msgId] = mkdir(mechDir,'solBF');
   %   if ~success
   %      error('Trouble making directory.');
   %   end
   %
   %   startCol = 1;
   %   for j = 1:numSolFiles
   %      endCol = min(2*dimBF,j*numBSolsPerFile);
   %      sol    = oldSol(startCol:endCol);
   %      solBFIdFile = ['solBF' filesep 'solId' ...
   %      sprintf(solFileIndexForm, startCol) '_' ...
   %      sprintf(solFileIndexForm,endCol)];
   %      save([mechDir filesep  solBFIdFile], 'sol');
   %      startCol = endCol+1;
   %   end
   %end


%   if exist('bfConstraint') == 1 && strcmp(bfConstraint,'adhLocation') == 1
%      %We also need to separate the finite elements that are inside and outside
%      % the adhesion sites. 
%      if strcmp(adhGeomType,'node') == 1
%         inMcfDOF = [];
%         for k = 1:length(indDomDOF)
%            if adhLabel(ceil(fs.mesh.p(2,indDomDOF(k))), ...
%               ceil(fs.mesh.p(1,indDomDOF(k)))) == 0
%               inMcfDOF = [inMcfDOF k];
%            end
%         end
%      elseif strcmp(adhGeomType,'convexHull') == 1
%         %If the adhesion geometry is given as the convex hull, we first build a 
%         % 'fem' structure for the adhesion geometry.
%         adhFem.geom  = adhGeom;
%         adhFem.mesh  = meshinit(adhFem);
%         adhFem.xmesh = meshextend(adhFem);
%
%         %Get the indices of the DOFS whose shape function are outside the adhesion
%         % sites and are also disconnected from the boundary.
%         [is,inMcfDOF] = postinterp(adhFem,fs.mesh.p(:,indDomDOF));
%      end
%
%      inAdhDOF = [1:length(indDomDOF)];
%      inAdhDOF(inMcfDOF) = [];
%      mcfDOF = indDomDOF(inMcfDOF);
%      adhDOF = indDomDOF(inAdhDOF);
%      dimMCF = length(mcfDOF);
%      dimAdh = dimBF-dimMCF;
%   end
