%This script file set up the model including the geometry, the mesh, the 
% coefficient functions in the equation and the boundary condition.

fprintf(1,'Set up the model : \n');

localStartTime = cputime;

%Read the cell image.
%cellImg = imread(firstImgFileName);

%Load model parameters
run([mechDir filesep 'modelPar']);
save([mechDir filesep 'lastSavedModParam.mat'],'modParam');

%rawDispFieldDir = [mechDir filesep 'rawDispField'];
%iDispFieldDir   = [reslDir filesep 'iDispField'];

%Set up the fem modle if it has not been done before.
if strcmp(bfFwdOpComputed,'none') == 1
   [modelFileList,modelFileImgIndex] = getFemModelFile(femModelDir);
   if isempty(modelFileList)
      fprintf(1,'Field boundary has not been drawn and saved yet. Run drawFieldBound first.\n');
      return;
   end

   answer = input(sprintf('Select Model fields (0 for all:%d-%d):',1,length(modelFileList)));
   if isempty(answer)
      selFields = 1:length(modelFileList);
   elseif length(answer) == 1 && answer == 0
      selFields = 1:length(modelFileList);
   else
      selFields = answer;
   end

   procStr = '';
   for ii = 1:length(selFields)
      jj = selFields(ii);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %Specify elasticity. We define the elasticity (Young's modulus)
      % to be proportional to the image intensity if inhomogeneous elasticity 
      % is chosen.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if exist('YModulVariation','var') && strcmp(YModulVariation,'inhomo')
         %First, get the image. We use the stacked image in the time average
         % window.
         startFrmNo = imgIndexOfDTimePts(jj)-firstImgIndex+1;
         if length(numAvgFrames) == numDTimePts
            curNumAvgFrames = numAvgFrames(jj);
         else
            curNumAvgFrames = numAvgFrames;
         end

         %Get the stacked images of the actin channel, pressumably the first channel.
         actImg = double(imread(imgFileList{1}{startFrmNo}));
         for k = 1:curNumAvgFrames-1
            actImg = actImg + double(imread(imgFileList{1}{k+startFrmNo}));
         end
         actImg = actImg/curNumAvgFrames;


         maxImgI = max(actImg(:));
         minImgI = min(actImg(:));
         if minImgI < maxImgI
            actImg = actImg/maxImgI;

            %Cut off a bright spot of very high intensity.
            fn.YModul = 'elImgParFun';
            fp.YModul = {{'x' 'y'} {actImg 1.5 2 0.1 'Gaussian' 5}};
         end
      end

      %for kk = 1:length(procStr)
      %   fprintf(1,'\b');
      %end
      procStr = sprintf('   Field %d ... ',jj);
      fprintf(1,procStr);

      clear femModel;

      s = load(modelFileList{jj});
      femModel = s.femModel;
      fieldBnd = femModel.fieldBnd;

      %Create the geometry for recovery from 'fieldBnd'. Since we need
      % fine mesh for the myosin contraction region bounded by 'bfDomPG'. So, 
      % we use 'rectmesh' and then define geometry with mesh as geometry.
      if strcmp(meshType,'rectMesh') && exist('bfConstraint') == 0
         curvL = [fieldBnd.x(fieldBnd.VI(1):fieldBnd.VI(2)) fieldBnd.y(fieldBnd.VI(1):fieldBnd.VI(2))].';
         curvT = [fieldBnd.x(fieldBnd.VI(2):fieldBnd.VI(3)) fieldBnd.y(fieldBnd.VI(2):fieldBnd.VI(3))].';
         curvR = [fieldBnd.x(fieldBnd.VI(3):fieldBnd.VI(4)) fieldBnd.y(fieldBnd.VI(3):fieldBnd.VI(4))].';
         curvB = [fieldBnd.x(fieldBnd.VI(4):end) fieldBnd.y(fieldBnd.VI(4):end)].';
         %msh   = rectmesh(curvL,curvB,curvT,curvR,50,[0:0.025:0.5 0.6:0.1:1]);
         msh   = rectmesh(curvL,curvT,curvR,curvB,recHin,recVin);

         %vertEI   = find(msh.e(3,:)==0); 
         %numEdges = length(vertEI);
         numEdges = max(msh.e(5,:));
         numSubDoms = 1;
         geom     = [];
      end

      %Separate the boundary segments into internal border and real boundaryies.
      borderE   = find(msh.e(6,:)~=0&msh.e(7,:)~=0);
      borderSeg = unique(msh.e(5,borderE));
      exBndE    = find(msh.e(6,:)==0|msh.e(7,:)==0);
      exBndSeg  = unique(msh.e(5,exBndE));
      numEdges  = length(exBndSeg);

      %Group subdomains and boundaries.
      ind    = ones(1,numSubDoms);
      bndInd = {};
      for k = 1:numEdges
         bndInd{k} = exBndSeg(k);
      end

      %Specify the boundary traction force.
      fn.BndDispx  = cell(1,numEdges);
      fn.BndDispy  = cell(1,numEdges);
      fp.BndDispx  = cell(1,numEdges);
      fp.BndDispy  = cell(1,numEdges);
      fn.BndTracFx = cell(1,numEdges);
      fn.BndTracFy = cell(1,numEdges);
      fp.BndTracFx = cell(1,numEdges);
      fp.BndTracFy = cell(1,numEdges);

      for k = 1:numEdges
         %Specify zero boundary displacement initially.
         fn.BndDispx{k} = 0;
         fn.BndDispy{k} = 0;

         %Set up link to boundary traction force function.
         fn.BndTracFx{k} = 'spBndTF';
         fn.BndTracFy{k} = 'spBndTF';
         fp.BndTracFx{k} = {{'s'} {0}};
         fp.BndTracFy{k} = {{'s'} {0}};

         %Specify the type of boundary conditions.
         BCTypes{k} = 'Dirichlet';
      end

      if ~isempty(borderE)
         bndInd{numEdges+1} =  borderSeg;

         fn.BndDispx{numEdges+1}  = [];
         fn.BndDispy{numEdges+1}  = [];
         fp.BndDispx{numEdges+1}  = [];
         fp.BndDispy{numEdges+1}  = [];

         fn.BndTracFx{numEdges+1} = 0;
         fn.BndTracFy{numEdges+1} = 0;
         fp.BndTracFx{numEdges+1} = [];
         fp.BndTracFy{numEdges+1} = [];

         BCTypes{numEdges+1} = 'Neumann';
      end

      options = elOptionsSet('EPType','YModulPRatio','BCType', BCTypes);

      %Specify the body force.
      fn.BodyFx = 'femBodyF';
      fn.BodyFy = 'femBodyF';
      fp.BodyFx = {{'x' 'y'} {[] 0}};
      fp.BodyFy = {{'x' 'y'} {[] 0}};

      if exist('fem','var') == 1
         clear fem; %Clear the 'fem' from simulation or the old 'fem'.
      end

      if isempty(geom)
         fem = elModelAssemble([],msh,options,fn,fp,ind,bndInd);
      else
         if exist('domMeshPar') == 1
            fem = elModelAssemble(geom,domMeshPar,options,fn,fp,ind,bndInd);
         else
            fem = elModelAssemble(geom,[],options,fn,fp,ind,bndInd);
         end
         msh = fem.mesh;
      end

      %Create the function space for domain force.
      %Name of any function in the space.
      %fs.dim = {'v1' 'v2'};
      fs.fem.dim = {'v'};

      fs.fem.shape = 1;
      fs.fem.bnd.h = 1;
      fs.fem.bnd.r = 0;

      fs.fem.mesh  = rectmesh(curvL,curvT,curvR,curvB,bfDomHin,bfDomVin);
      fs.fem.geom  = fs.fem.mesh;

      fs.fem.equ.ind = ones(1,numSubDoms);
      fs.fem.xmesh = meshextend(fs.fem);

      %Find the indices of the DOFs whose shape functions have support disconnected 
      % from the boundary.
      N = assemble(fs.fem,'out','N');
      indDomDOF = find(full(sum(N,1))==0);

      %The Degree of Freedom vector for the basis or shape function in the finite
      % element space.
      %Dimension of the function space where the basis can be nonzero on the 
      % boundary.
      dimFS  = flngdof(fs.fem);

      %Dimention of the function space where the basis function is kept zero on the
      % boundary.
      dimBF = length(indDomDOF);

      fs.dimFS     = dimFS;
      fs.dimBF     = dimBF;
      fs.indDomDOF = indDomDOF;

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %                  Four boundary edges
      % The geometry and mesh, 'msh' created by 'rectmesh' has four boundary edges
      % parameterized by arclength. So, first identify the four edges.
      % We use a structure named 'edge' to store all these mesh structure info and
      % displacement. It has the following fields:
      % 'vertEI' : Index into 'msh.e(3,:)' whose arclength is 0. It identifies the
      %            vertex of the edge.
      % 'I'      : Index into 'msh.p' to get the real coordinates of boundary points.
      % 'bndEI'  : The index of boundary elements that belong to one edge.
      % 'endI    : The index of the boundary element that is at the end of one
      %            edge.
      % 'arcLen' : The arclength of each edge.
      % 'bndP'   : The coordinates of boundary points on each edge.
      % 'bndS'   : The arclength parameters of the boundary points on each edge.
      % 'ppX'    : The pp-form of the spline interpolation of X-coordinate of 'bndP'.
      % 'ppY'    : The pp-form of the spline interpolation of Y-coordinate of 'bndP'.
      % 'dispV'  : The coordinates of the base and the point end of the displacement
      %            vectors on each edge in the formate [x0 y0 x1 y1].
      % 'U1'     : The first and 
      % 'U2'     : the second components of the displacement vectors on each edge.
      % 'UC1'    : For debugging. Has the same structure as 'edgeU1'.
      % 'UC2'    : For debugging. Has the same structure as 'edgeU2'.
      % 'ppU1'   : The pp-form of the spline interpolation of 'U1'.
      % 'ppU2'   : The pp-form of the spline interpolation of 'U2'.
      for k = 1:numEdges
         %Extract the arclenth parameters.
         edge(k).bndEI = find(msh.e(5,:)==k);
         edge(k).bndS  = msh.e(3,edge(k).bndEI);

         [edge(k).bndS,sortI] = sort(edge(k).bndS);

         edge(k).endI   = edge(k).bndEI(sortI(end));
         edge(k).bndS   = [edge(k).bndS msh.e(4,edge(k).endI)];
         edge(k).arcLen = edge(k).bndS(end);

         %Get the index of the boundary points.
         edge(k).I = msh.e(1,edge(k).bndEI);
         edge(k).I = [edge(k).I(sortI) msh.e(2,edge(k).endI)];

         edge(k).bndP  = [msh.p(1,edge(k).I); msh.p(2,edge(k).I)];
         bndSKnt = augknt(edge(k).bndS,2);
         edge(k).ppX = spapi(bndSKnt,edge(k).bndS, edge(k).bndP(1,:));
         edge(k).ppY = spapi(bndSKnt,edge(k).bndS, edge(k).bndP(2,:));
         %edge(k).ppX = spline(edge(k).bndS, edge(k).bndP(1,:));
         %edge(k).ppY = spline(edge(k).bndS, edge(k).bndP(2,:));
      end

      femModel.fem      = fem;
      femModel.fs       = fs;
      femModel.options  = options;
      femModel.fn       = fn;
      femModel.fp       = fp;
      femModel.curvL    = curvL;
      femModel.curvT    = curvT;
      femModel.curvR    = curvR;
      femModel.curvB    = curvB;
      femModel.BCTypes  = BCTypes;
      femModel.geom     = geom;
      femModel.ind      = ind;
      femModel.bndInd   = bndInd;
      femModel.numEdges = numEdges;
      femModel.edge     = edge;
      femModel.edgeMsh   = msh; %Mesh used for identify boundary edgies.

      femModelFile = [femModelDir filesep 'femModel' ...
         sprintf(imgIndexForm,modelFileImgIndex(jj)) '.mat'];
      save(femModelFile,'femModel');

      %save([mechDir filesep 'femId'],'fem');
      %load([mechDir filesep 'femId'],'fem');
      %load(femModelFile);
      %msh = fem.mesh;
      fprintf(1,'Done in %5.3f sec.\n',cputime-localStartTime);
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get preprocessed experimental data such as calculating boundary displacement. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(1,'Identify domain data points and calculate boundary displacement ... \n');
answer = input('Select time steps (0 for all):');
if isempty(answer) | answer == 0
   selTimeSteps = 1:numDTimePts;
else
   selTimeSteps = answer;
end

for ii = 1:length(selTimeSteps)
   jj = selTimeSteps(ii);

   localStartTime = cputime;
   fprintf(1,'   Time Step: %d ... ',jj);

   imgIndex = imgIndexOfDTimePts(jj);
   % Load the raw field
   rawDispFieldFileName = ['rawDispField' sprintf(imgIndexForm,imgIndex) '.mat'];
   rawDispFieldFile     = [rawDispFieldDir filesep rawDispFieldFileName];

   if exist(rawDispFieldFile,'file')
      s = load(rawDispFieldFile);
      rawDispField = s.rawDispField;
   else
      fprintf(1,'Raw vector field has not been loaded and saved yet. Run loadRawField first.\n');
      return;
   end

   %Load the field boundary.
   if strcmp(isFieldBndFixed,'yes')
      femModelFile = [femModelDir filesep 'femModel' sprintf(imgIndexForm,0) '.mat'];
   else
      femModelFile = getFemModelFile(femModelDir,imgIndex,imgIndexForm);
   end

   clear femModel;

   if exist(femModelFile,'file')
      s = load(femModelFile);
      femModel = s.femModel;
   else
      fprintf(1,'Fem model file is missing.\n');
      return;
   end

   fieldBnd = femModel.fieldBnd;
   fem      = femModel.fem;
   numEdges = femModel.numEdges;
   edge     = femModel.edge;

   %Load the interpolated displacement field.
   iDispFieldFileName = ['iDispField' sprintf(imgIndexForm,imgIndex) '.mat'];
   iDispFieldFile     = [iDispFieldDir filesep iDispFieldFileName];
   s = load(iDispFieldFile);
   iDispField = s.iDispField;

   [is,pe] = postinterp(fem,iDispField.p.');
   iDispField.iOutMesh    = pe; % 'pe': index of points outside 'msh'.
   iDispField.iInMesh     = 1:size(iDispField.p,1);
   iDispField.iInMesh(pe) = [];
   iDispField.numDP       = length(iDispField.iInMesh);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %Calculate the displacements of the boundary edges given by 'msh'. See
   % 'help meshinit' for information about the MESH structure in FEMLAB.
   if strcmp(trackMethod,'simu') && ~strcmp(dataToUse,'smooth')
      smoothFactor = input('Select smooth factor (0-5): ');
      if isempty(smoothFactor)
         smoothFactor = 0;
      end

      edgD = rawDispField.edgD;
      for kk = 1:numEdges
         edgS = linspace(edgD(kk).s(1),edgD(kk).s(end), ...
            ceil(length(edgD(kk).s)/(smoothFactor+2)));
         sKnot = augknt(edgS,4);
         %sKnot = augknt(edgD(kk).s(1:2:end),4);
         edgD(kk).ppU1 = spap2(sKnot,4,edgD(kk).s,edgD(kk).U1.');
         edgD(kk).ppU2 = spap2(sKnot,4,edgD(kk).s,edgD(kk).U2.');
         edgD(kk).smoothFactor = smoothFactor;
      end
   else
      edgCorLen = iDispField.edgCorLen;
      for k = 1:numEdges
         rawDispV = [rawDispField.p(:,2:-1:1) rawDispField.p(:,2:-1:1)+rawDispField.v(:,2:-1:1)];
         edgD(k).dispV  = vectorFieldSparseInterp(rawDispV, ...
            edge(k).bndP(2:-1:1,:).',2*edgCorLen,edgCorLen,[]);
         %edgeU{k}  = vectorFieldInterp(rawDispV{jj}, ...
         %   edgeP{k}(2:-1:1,:).',edgCorLen,[]);
         edgD(k).U1 = edgD(k).dispV(:,4) - edgD(k).dispV(:,2);
         edgD(k).U2 = edgD(k).dispV(:,3) - edgD(k).dispV(:,1);

         %For debugging:
         %edgeU1{jj,k} = edgeU1{jj,k}/2; 
         %edgeU2{jj,k} = edgeU2{jj,k}/2;

         %Create spline interpolation of the edge displacement using arclength
         % parameter stored in 'msh.e(3:4,:)'.
         numInd = find(~isnan(edgD(k).U1));
         edgD(k).s = edge(k).bndS(numInd);
         bndSKnt = augknt(edge(k).bndS(numInd),2);
         edgD(k).ppU1 = spapi(bndSKnt,edge(k).bndS(numInd), edgD(k).U1(numInd).');
         %edgD(k).ppU1 = spline(edge(k).bndS(numInd), edgD(k).U1(numInd).');
         numInd = find(~isnan(edgD(k).U2));
         bndSKnt = augknt(edge(k).bndS(numInd),2);
         edgD(k).ppU2 = spapi(bndSKnt,edge(k).bndS(numInd), edgD(k).U2(numInd).');
         %edgD(k).ppU2 = spline(edge(k).bndS(numInd), edgD(k).U2(numInd).');
      end
   end

   %Save the edge mesh information and displacements data.
   iDispField.edgD = edgD;
   save(iDispFieldFile,'iDispField');

   fprintf(1,'Done in %5.3f sec.\n', cputime-localStartTime);
end

%Save the displacements on the boundary.
%save([reslDir filesep  'edgeDisp'],'numEdges','edgePPx','edgePPy', ...
%   'edgeS','edgeP','edgeU1','edgeU2');

%Save the displacements on the data points.
%save([reslDir filesep  'dataDisp'],'dataPx','dataPy','dataU1','dataU2', ...
%   'rawDataP');



%Some old code.
%   if strcmp(meshType,'rectMesh') && exist('bfConstraint') == 0
%      curvL = [fieldBnd.x(fieldBnd.VI(1):fieldBnd.VI(2)) fieldBnd.y(fieldBnd.VI(1):fieldBnd.VI(2))].';
%      curvT = [fieldBnd.x(fieldBnd.VI(2):fieldBnd.VI(3)) fieldBnd.y(fieldBnd.VI(2):fieldBnd.VI(3))].';
%      curvR = [fieldBnd.x(fieldBnd.VI(3):fieldBnd.VI(4)) fieldBnd.y(fieldBnd.VI(3):fieldBnd.VI(4))].';
%      curvB = [fieldBnd.x(fieldBnd.VI(4):end) fieldBnd.y(fieldBnd.VI(4):end)].';
%      %msh   = rectmesh(curvL,curvB,curvT,curvR,50,[0:0.025:0.5 0.6:0.1:1]);
%      msh   = rectmesh(curvL,curvT,curvR,curvB,recHin,recVin);
%
%      %vertEI   = find(msh.e(3,:)==0); 
%      %numEdges = length(vertEI);
%      numEdges = max(msh.e(5,:));
%      numSubDoms = 1;
%      geom     = [];
%   else
%      curve = geomspline([fieldBnd.x fieldBnd.y].', ...
%      'splinemethod','centripetal','closed','auto');
%      geom = geomcoerce('solid',{curve});
%      numSubDoms = geominfo(geom,'Out',{'nmr'});
%
%      if exist('bfConstraint') == 1 && strcmp(bfConstraint,'adhLocation') == 1
%         if strcmp(isAdhSegmented,'no') == 1
%            fprintf(1,'   Get the geometry of the adhesion site.\n');
%            %If there is the constraint of adhesion locations, identify the 
%            % adhesion sites using the cluster and segmentation technique.
%            adhSeg = imClusterSeg(double(adhImg),1,'method','kmeans', ...
%            'k_cluster',numClusters);
%            adhSeg(find(adhSeg<numClusters)) = 0;
%
%            %We only consider adhesion sites that are inside the geometry domain.
%            % To get rid of the adhesion sites that are outside, we first create 
%            % the geometry and meshing of the domain.
%            clear fem;
%            fem.geom = geom;
%            if exist('domMeshPar') == 1
%               fem.mesh = meshinit(fem,domMeshPar{:});
%            else
%               fem.mesh = meshinit(fem);
%               fem.mesh = meshrefine(fem);
%            end
%
%            fem.xmesh = meshextend(fem);
%
%            %To get rid of the adhesion pixels that are outside the domain of
%            % interest, first create the pixel list of the whole image.
%            %Number of pixels in the horizontal direction
%            numPixelsH = size(adhSeg,2); 
%            %Number of pixels in the vertical direction
%            numPixelsV = size(adhSeg,1); 
%            [pixelH,pixelV] = meshgrid([1:numPixelsH],[1:numPixelsV]);
%            pixelH = reshape(pixelH,1,length(pixelH(:)));
%            pixelV = reshape(pixelV,1,length(pixelV(:)));
%
%            [is,pe] = postinterp(fem,[pixelH;pixelV]);
%            adhSeg(pe) = 0;
%
%            %Get rid of small adhesion sites. 
%            adhLabel   = bwlabel(adhSeg,8); 
%            adhStat    = regionprops(adhLabel,'Area');
%            adhArea    = [adhStat.Area];
%            smAdhLabel = find(adhArea<smAdhThreshold);
%            for k = 1:length(smAdhLabel)
%               adhSeg(find(adhLabel==smAdhLabel(k))) = 0;
%            end
%
%            %Relabel the adhesion segmentation.
%            adhLabel = bwlabel(adhSeg,8); 
%
%            save([mechDir filesep 'adhSeg'],'adhLabel','adhSeg');
%         else
%            fprintf(1,'  Load the segmented adhesion sites.\n');
%            load([mechDir filesep 'adhSeg']);
%         end
%
%         if strcmp(isAdhGeomBuilt,'no') == 1
%            %We want to have fine enough mesh at the adhesion site. We propose
%            % two ways to achieve this. The first approach is to select a set of
%            % seeding points on the adhesion sites and force them to be the nodes
%            % of meshing. The second approach is to find the convex hull of the
%            % adhesion site and use the 'hmax' property of the 'meshinit' command 
%            % to assign different maximun meshing size on the subdomain of 
%            % adhesion sites.
%
%            if strcmp(adhGeomType,'node') == 1
%               %Get the set of pixels for each adhesion site
%               adhPixelList = {};
%               for k = 1:max(adhLabel(:))
%                  [adhPixelList{k}.y,adhPixelList{k}.x] = find(adhLabel==k);
%
%                  %Get the boundary pixels for the singled out adhesion site.
%                  %adhBnd{k} = bwperim(adhSeg,4);
%                  %[bndy,bndx] = find(adhBnd{k}~=0);
%                  %adhBnd{k} = sortEdgePixels([bndx bndy]);
%                  %adhSeg(:) = 0;
%               end
%
%               adhNodes = {};
%               numNodes = 0;
%               for k = 1:length(adhPixelList)
%                  %Get the x coordinate of the leftmost and the rightmost pixel.
%                  xmin = min(adhPixelList{k}.x);
%                  xmax = max(adhPixelList{k}.x);
%
%                  %Calculate the number of vertical grid lines in the x 
%                  % (horizontal) direction.
%                  numVGridLines = floor((xmax-xmin)/adhMeshDx+0.5)+1;
%
%                  %Calculate the x coordinates of these vertical grid lines
%                  if numVGridLines < 3 
%                     vGridLinesX = ceil((xmin+xmax)/2);
%                     numVGridLines = 1;
%                  else
%                     vGridLinesX = ceil(linspace(xmin,xmax,numVGridLines));
%                  end
%
%                  %Along each line, we choose evenly distributed pixels as the seed
%                  % points for meshing.
%                  for jj = 1:numVGridLines
%                     %Get all the adhesion pixels along each line
%                     ind = find(adhPixelList{k}.x == vGridLinesX(jj));
%
%                     %Get the top and bottom pixel on the line.
%                     ymin = min(adhPixelList{k}.y(ind));
%                     ymax = max(adhPixelList{k}.y(ind));
%
%                     %Calculate the number of seed points along this line.
%                     numPtsY = floor((ymax-ymin)/adhMeshDy+0.5)+1;
%                     if numPtsY < 3 
%                        pixelY  = (ymin+ymax)/2;
%                        numPtsY = 1;
%                     else
%                        pixelY = linspace(ymin,ymax,numPtsY);
%                     end
%
%                     %Add the seed points to the set of nodes.
%                     for nl = 1:length(pixelY)
%                        numNodes = numNodes+1;
%                        adhNodes{numNodes} = point2(vGridLinesX(jj),pixelY(nl));
%                     end
%                  end
%               end
%
%               save([mechDir filesep 'adhGeom'],'adhNodes');
%               %%Mesh the whole domain.
%               %clear fem;
%               %fem.geom = geom;
%               %if exist('domMeshPar') == 1
%               %   fem.mesh = meshinit(fem,domMeshPar{:});
%               %else
%               %   fem.mesh = meshinit(fem);
%               %   fem.mesh = meshrefine(fem);
%               %end
%
%               %%Get rid of the adhesion meshing points that are outside the domain of
%               %% interest.
%               %fem.xmesh = meshextend(fem);
%               %[is,pe] = postinterp(fem,adhMsh.p);
%               %adhMsh.p(:,pe) = [];
%
%               %%Get rid of meshing points that are on the adhesion sites.
%               %msh = fem.mesh;
%               %ind = [];
%               %for k = 1:size(msh.p,2)
%               %   if adhLabel(ceil(msh.p(2,k)),ceil(msh.p(1,k))) ~= 0
%               %      ind = [ind k];
%               %   end
%               %end
%               %msh.p(:,ind) = [];
%
%               %Combine 'msh' and 'adhMsh'.
%               %msh2.p = [msh.p adhMsh.p];
%               %msh    = meshenrich(msh2);
%               %else
%               %   numEdges = geominfo(geom,'Out',{'nbs'});
%            elseif strcmp(adhGeomType,'convexHull') == 1
%               %First, find the convex hull of each adhesion.
%               adhStat         = regionprops(adhLabel,'ConvexHull');
%               fineAdhConvHull = {adhStat.ConvexHull};
%
%               %Generate the boundary curve of each adhesion site from 
%               % 'finAdhConvHull'.
%               adhBndCurve = {};
%               adhConvHull = {};
%               for k = 1:length(fineAdhConvHull)
%                  %First Eliminate some vertices on the 'adhConvHull' so that the 
%                  % meshes generated later will not be too small.
%                  adhConvHull{k} = fineAdhConvHull{k}(1,:);
%                  arcLen  = [0; sqrt(diff(fineAdhConvHull{k}(:,1)).^2 + ...
%                  diff(fineAdhConvHull{k}(:,2)).^2)];
%                  for jj = 1:length(arcLen)-1
%                     arcLen(jj+1) = arcLen(jj+1)+arcLen(jj);
%                  end
%
%                  smArcLenThreshold = min(3*sqrt(smAdhThreshold),arcLen(end)/4);
%                  bigArcLen         = smArcLenThreshold;
%                  for jj = 1:length(arcLen)-1
%                     if arcLen(jj+1) > bigArcLen
%                        adhConvHull{k} = [adhConvHull{k}; ...
%                        (fineAdhConvHull{k}(jj+1,:)*(bigArcLen-arcLen(jj))+...
%                        fineAdhConvHull{k}(jj,:)*(arcLen(jj+1)-bigArcLen))/...
%                        (arcLen(jj+1)-arcLen(jj))];
%                        bigArcLen = bigArcLen + smArcLenThreshold;
%                     end
%                  end
%                  if norm(adhConvHull{k}(end,:)-fineAdhConvHull{k}(end,:)) < ...
%                     smArcLenThreshold/2
%                     adhConvHull{k}(end,:) = fineAdhConvHull{k}(end,:);
%                  else
%                     adhConvHull{k}(end+1,:) = fineAdhConvHull{k}(end,:);
%                  end
%
%                  adhBndCurve{k} = geomspline(adhConvHull{k}.', ...
%                  'splinemethod','centripetal','closed','auto');
%               end
%
%               adhGeom = geomcoerce('solid',adhBndCurve);
%               save([mechDir filesep 'adhGeom'],'adhGeom','adhBndCurve');
%            end
%         else
%            fprintf(1,'  Load the geometry of the adhesion sites.\n');
%            load([mechDir filesep 'adhGeom']);
%         end
%
%         if strcmp(adhGeomType,'node') == 1
%            geom = geomcoerce('solid',{curve adhNodes{:}});
%         else
%            geom = geom+adhGeom;
%         end
%
%         %Get the number of subdomains.
%         numSubDoms = geominfo(geom,'Out',{'nmr'});
%
%         %Specify a smaller enought meshing size for adhesion sites. Note:
%         % Subdomain No. 1 should be the main domain surrounding the adhesion
%         % sites.
%         subHmax = zeros(2,numSubDoms-1);
%         for k = 2:numSubDoms
%            subHmax(:,k-1) = [k;sqrt(smAdhThreshold)];
%         end
%
%         %Get the global 'hmax' if it is specified in 'domMeshPar'.
%         %The position of the 'hmax' property.
%         if exist('domMeshPar')
%            hmaxId = find(strcmp(domMeshPar,'hmax'));
%         else
%            hmaxId = [];
%         end
%
%         if isempty(hmaxId)
%            %A gloabl 'hmax' is not specified in 'domMeshPar'.
%            domMeshPar{end+1} = 'hmax';
%            domMeshPar{end+1} = {sqrt(smAdhThreshold) [] [] subHmax};
%         else
%            %Merge the global 'hmax' with the adhesion 'hmax'.
%            domMeshPar{hmaxId+1} = {domMeshPar{hmaxId+1} [] [] subHmax};
%         end
%      end
%
%      %Mesh the geometry
%      if exist('domMeshPar')
%         msh = meshinit(geom,domMeshPar{:});
%      else
%         fem.geom = geom;
%         fem.mesh = meshinit(fem);
%         msh = meshrefine(fem);
%      end
%      geom     = [];
%   end
