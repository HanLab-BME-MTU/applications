function reconstructMovieIntracelluarForce(movieData,varargin)
% reconstructMovieIntracelluarForce reconstructs intracellular forces out
% of speckle flow. Adaped from Ji, Lin's code in contMechModel2. This
% function this time tries to use Matlab's PDE toolbox instead of using
% comsol 3.5 functions, which changed after version 4.0 significantly
% Sangyoon Han Jan 2015
%% ----------- Input ----------- %%

%Check input
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('movieData', @(x) isa(x,'MovieData'));
ip.addOptional('paramsIn',[], @isstruct);
ip.addParameter('saveOnlyMesh',false,@islogical);
ip.parse(movieData,varargin{:});
saveOnlyMesh=ip.Results.saveOnlyMesh;

%Get the indices of any previous flow tracking processes                                                                     
iFlowProc = movieData.getProcessIndex('FlowTrackingProcess',1,0);
p = movieData.getProcess(iFlowProc).funParams_;

% iProc = movieData.getProcessIndex('FlowAnalysisProcess',1,0);
% %If the process doesn't exist, create it
% if isempty(iProc)
%     iProc = numel(movieData.processes_)+1;
%     movieData.addProcess(FlowAnalysisProcess(movieData,...
%         movieData.outputDirectory_));                                                                                                 
% end
% 
% flowAnProc = movieData.processes_{iProc};
% %Parse input, store in parameter structure
% p = parseProcessParams(flowAnProc,paramsIn);
% if strcmp(p.FlowProcess,'SpeckleTrackingProcess')
%     p.FlowProcess = 'FlowTrackingProcess';
% end

%% --------------- Initialization ---------------%%
if feature('ShowFigureWindows')
    wtBar = waitbar(0,'Initializing...','Name','Intracellular force reconstruction ...');
end

% Reading various constants
nFrames = movieData.nFrames_;

% Test the presence and output validity of the speckle tracking process
iSpecDetProc =movieData.getProcessIndex('SpeckleDetectionProcess',1,1);     
if isempty(iSpecDetProc)
    error(['Speckle detection has not yet been performed'...
    'on this movie! Please run first!!']);
end        
%Check that there is a valid output
specDetProc = movieData.processes_{iSpecDetProc};
% if ~specDetProc.checkChannelOutput(p.ChannelIndex)
%     error(['Each channel must have speckles! ' ...
%         'Please apply speckle detection to all needed channels before '...
%         'running flow analysis!'])
% end
p.MaskChannelIndex = specDetProc.funParams_.MaskChannelIndex;

% Create mask directory if several masks need to be merged
if length(p.MaskChannelIndex) >1
    %Get the indices of any previous mask intersection process
    iMaskProc = movieData.getProcessIndex('MaskIntersectionProcess',1,0);
else  
    iMaskProc =movieData.getProcessIndex('MaskRefinementProcess',1,1);
end

if isempty(iMaskProc)
    error(['Mask refinement has not yet been performed '...
        'on this movie! Please run first!!']);
end
%Check that there is a valid output
maskProc = movieData.processes_{iMaskProc};
% if ~maskProc.checkChannelOutput(p.ChannelIndex)
%     error(['Each channel must have masks !' ...
%         'Please apply mask refinement to all needed channels before'...
%         'running speckle tracking!'])
% end

% Test the presence and output validity of the speckle tracking process
% iFlowProc =movieData.getProcessIndex(p.FlowProcess,1,1);     
% if isempty(iFlowProc)
%     error([eval([p.FlowProcess '.getName']) ' has not yet been performed'...
%     ' on this movie! Please run first!!']);
% end        

%Check that there is a valid output
% flowProc = movieData.processes_{iFlowProc};
flowProc = movieData.processes_{iFlowProc};
% if ~flowProc.checkChannelOutput(p.ChannelIndex)
%     error(['Each channel must have flow! Please apply '...
%         eval([p.FlowProcess '.getName']) ' to all needed channels before '...
%         'running flow analysis!'])
% end
% this input part should be updated once the algorithm is completed.
% imgIndexForm  = sprintf('%%.%dd',length(num2str(no)));
%% Output folder setup
outputDir = [movieData.outputDirectory_ filesep 'IFMPackage'];
mkClrDir(outputDir)
femModelDir = [outputDir filesep 'femModel'];
mkClrDir(femModelDir)
iDispFieldDir = [outputDir filesep 'iDispField'];
mkClrDir(iDispFieldDir)
meshTifPath = [outputDir filesep 'meshTif'];
mkClrDir(meshTifPath)
%% Import the flow data per each frame - the flow should've been calculated
% between only the two adjacent frames
%% --------------- Kinetic analysi ---------------%%% 

disp('Starting reconstructing forces...')
% Format string for zero-padding file names
% Set up the input directories
nChan=numel(movieData.channels_);
outputDir=cell(1,nChan);
fString = ['%0' num2str(floor(log10(nFrames))+1) '.f'];
numStr = @(frame) num2str(frame,fString);
% Anonymous functions for reading input/output
logMsg = @(chan) ['Please wait, analyzing flow for channel ' num2str(chan)];
outFile=@(chan,frame) [outputDir{chan} filesep 'intraForces_' numStr(frame) '.mat'];

numSubDoms = 1; % this should be designated in the dialog box.
for i=1:numel(p.ChannelIndex)
    iChan = p.ChannelIndex(i);
    % Log display
%     channelLog{i} = sprintf('Channel %g: %s\n',iChan,inFilePaths{1,iChan});
    disp(logMsg(iChan))
    
    if ishandle(wtBar), waitbar(0,wtBar,'Loading masks and images...'); end
%     maskNames = maskProc.getOutMaskFileNames(iChan);
%     inMask=@(frame) [flowAnProc.inFilePaths_{2,iChan}...
%         filesep maskNames{1}{frame}];
    mask=true([movieData.imSize_ nFrames]);
    stack=zeros([movieData.imSize_ nFrames]);
    for jj = 1:nFrames
%         mask(:,:,jj) = logical(imread(inMask(jj)));
        mask(:,:,jj) = maskProc.loadChannelOutput(iChan,jj);
        stack(:,:,jj) = movieData.channels_(iChan).loadImage(jj);
    end

    % Load candidates and generate Nx3 matrices with position and intensity
    % Replace fsmTrackFillSpeckleList
    if ishandle(wtBar), waitbar(0,wtBar,['Loading flow for channel ' num2str(iChan)']); end
    
    flow = flowProc.loadChannelOutput(iChan,'output','flow');
    flow=flow(1:end-1);
    % Remove NaNs 
    for jj=find(~cellfun(@isempty,flow))
        flow{jj}=flow{jj}(~isnan(flow{jj}(:,3)),:);
    end           

    % Interpolate field
    if ishandle(wtBar), waitbar(.25,wtBar,['Interpolating flow for channel ' num2str(iChan)']); end
    [Md,~,~,~,~] =  analyzeFlow(flow,p.timeWindow+1,p.maxCorLength,...
        'noise',1,'error',1);
    
   %% Mesh generation (decsg and initmesh)
    % It might be good to use the mask from mask process
    % and determine free boundary vs. inner boundary
    % by looking at its straightness
    % Take the binary image into geometry with 4 edges (3 inner + 1 free)
    % mask is mask(:,:,j)
%     saveOnlyMesh = 1; %intermediate step...
    iiformat = ['%.' '3' 'd'];
    display('Meshing and calculating the displacements on boundary edges...')
    for jj = 1:nFrames-1
        % Adjust mask boundaries with speckle boundaries
        curFlow = flow{jj};
 
%         B{1}(ymaxIdx,1) = RightLowerCorner(2); % for bottom
%         B{1}(xmaxIdx,2) = RightLowerCorner(1); % for right
%         B{1}(yminIdx,1) = LeftUpperCorner(2); % for top
%         B{1}(xminIdx,2) = LeftUpperCorner(1); % for left
        %% geometry and mesh
        minImgSize = 5; % edge length should be more than 5 pixel.
        if saveOnlyMesh
            [msh,borderE,borderSeg,exBndE,exBndSeg,numEdges,bndInd,ind,hFig] = getMeshFromMask(movieData,jj,curFlow, mask(:,:,jj),minImgSize,numSubDoms,1);
            I = getframe(hFig);
            imwrite(I.cdata, strcat(meshTifPath,'/meshTif',num2str(jj,iiformat),'.tif'));
            close(hFig)
            continue
        else
            [msh,borderE,borderSeg,exBndE,exBndSeg,numEdges,bndInd,ind] = getMeshFromMask(movieData,jj,curFlow, mask(:,:,jj),minImgSize,numSubDoms);
        end            
        %% PDE definition for continuum mechanics, plane stress
        % Now I have a mesh information, with which I can solve for forward
        % solution for each basis function
        % I am using PDE structural mechanics, plane stress
        % governing equation:
        %           -div(c*grad(u)) = 0
        % where c11 = [2*G+mu 0; 0 G], c12 = [0 G; mu 0], c21=[0 G;mu 0],
        % c22 = [G 0; 0 2*G+mu];
        % Material properties
        E = 1e1; % Elastic modulus, Pa = N/m^2. I need to think about how to make it spatially heterogeneous later
        Nu = 0.29; % Poisson's ratio for actin network. need to find out more correct one
        G = E/(2*(1+Nu)); % Shear modulus, Pa
        mu = 2*G*Nu/(1-Nu); % mu
        c11 = [2*G+mu 0; 0 G];
        c12 = [0 G; mu 0];
        c21=[0 G;mu 0];
        c22 = [G 0; 0 2*G+mu];
        c2d = [c11 c12;c21 c22]; % 2Nx2N flattened coefficient c
        c = c2d(:); % vectorized c
        a = 0; %a=0 for structural mechanics
        f = [0 0]';
        % I have to do something like this: fem = elModelAssemble([],msh,options,fn,fp,ind,bndInd);
        % where fn and fp contains boundary conditions etc.. 

        %% Boundary condition definition
        [fn,fp,BCTypes,bndInd] = initializeCMBoundaryCondition(numEdges,borderE,borderSeg,bndInd);
        options = elOptionsSet('EPType','YModulPRatio','BCType', BCTypes);
%         fem = elModelAssemble([],msh,options,fn,fp,ind,bndInd);
        fem = elModelAssemblePDE([],msh,options,fn,fp,ind,bndInd);

        %% function space creation for domain force.
        fs = constructFunctionSpace(msh,numSubDoms);
        %% Four boundary edges
        % The geometry and mesh, 'msh' created by 'getMeshFromMask' has four boundary edges
        % parameterized by arclength. So, first identify the four edges from msh.e.
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
        edge = edgeFromMsh(msh,numEdges);

        % saving all these into one struct variable femModel
        femModel.fem      = fem;
        femModel.fs       = fs;
        femModel.options  = options;
        femModel.fn       = fn;
        femModel.fp       = fp;
%         femModel.curvL    = curvL;
%         femModel.curvT    = curvT;
%         femModel.curvR    = curvR;
%         femModel.curvB    = curvB;
        femModel.BCTypes  = BCTypes;
        femModel.geom     = msh;
        femModel.ind      = ind;
        femModel.bndInd   = bndInd;
        femModel.numEdges = numEdges;
        femModel.edge     = edge;
        femModel.edgeMsh   = msh; %Mesh used for identifying boundary edges.

        femModelFile = [femModelDir filesep 'femModel' num2str(jj) '.mat'];
        save(femModelFile,'femModel');

%         % Using speckle postions for mesh generation
%         dt=delaunayTriangulation(flow(:,2),flow(:,1));
%         triplot(dt)
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Get preprocessed experimental data such as calculating boundary displacement. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       localStartTime = cputime;
       fprintf(1,'   Time Step: %d ... ',jj);

%        imgIndex = imgIndexOfDTimePts(jj);
       % Load the raw field
       rawDispField = flow{jj};

       %Load the interpolated displacement field.
       iDispField = Md{jj};

%        [is,pe] = postinterp(fem,iDispField.p.');
%        [is,pe] = pdeintrp(msh.p,msh.t,iDispField.p.');
%        iDispField.iOutMesh    = pe; % 'pe': index of points outside 'msh'.
%        iDispField.iInMesh     = 1:size(iDispField.p,1);
%        iDispField.iInMesh(pe) = [];
%        iDispField.numDP       = length(iDispField.iInMesh); % I don't
%        think saving these parameters matters too much

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %Calculate the displacements of the boundary edges given by 'msh'. See
       % 'doc initmesh' for information about the MESH structure.
      edgCorLen = p.maxCorLength;
      for k = 1:numEdges
%          rawDispV = [rawDispField.p(:,2:-1:1) rawDispField.p(:,2:-1:1)+rawDispField.v(:,2:-1:1)];
         rawDispV = rawDispField; 
         edgD(k).dispV  = vectorFieldSparseInterp(rawDispV, ...
            edge(k).bndP(2:-1:1,:).',3*edgCorLen,edgCorLen,[]); % I need to check if this is more exact than interpolating vec.
         jj = 3;
         while sum(~isnan(edgD(k).dispV(:,4)))<2
             jj = jj+2;
             edgD(k).dispV  = vectorFieldSparseInterp(rawDispV, ...
                edge(k).bndP(2:-1:1,:).',jj*edgCorLen,edgCorLen,[]); % I need to check if this is more exact than interpolating vec.
         end
         edgD(k).U1 = edgD(k).dispV(:,4) - edgD(k).dispV(:,2);
         edgD(k).U2 = edgD(k).dispV(:,3) - edgD(k).dispV(:,1);

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

       %Save the edge mesh information and displacements data.
       iDispField_pos = iDispField(2:-1:1,:);
       iDispField_vec = iDispField(4:-1:3,:)-iDispField(2:-1:1,:);
       clear iDispField
       iDispField.p = iDispField_pos;
       iDispField.v = iDispField_vec;
       iDispField.edgD = edgD;
       iDispFieldFile = [iDispFieldDir filesep 'iDispField' num2str(jj)];
       save(iDispFieldFile,'iDispField');

       fprintf(1,'Done in %5.3f sec.\n', cputime-localStartTime);
    end

    %% calFwdOpBF - compatible section
    % calFwdOpBF computes the matrix approximation to the forward linear
    % operator for the body force.
    % but now we are calculating these using PDE toolbox instead of FEMLAB
    fprintf(1,'Calculating the forward operator for the Body Force :\n');
    %Step 1: Construct the matrix approximation to the forward operator.
    %We use each basis function as the body force to solve our
    % continuum mechanics system. The solution gives us each column of the matrix.
    dimFS     = length(msh.p);
%     dimBF     = fs.dimBF;
%     indDomDOF = fs.indDomDOF;
%     coefFS    = zeros(dimFS,1);
%     solFileIndexForm = sprintf('%%.%dd',length(num2str(2*dimBF)));
%     if rem(2*dimBF,numBSolsPerFile) == 0
%      numSolFiles = 2*dimBF/numBSolsPerFile;
%     else
%      numSolFiles = ceil(2*dimBF/numBSolsPerFile);
%     end

     %'k' is the index of 'coefFS' whose corresponding basis function is zero
     % on the boundary.
%      k = indDomDOF(j);
     for k=1:dimFS
         coefFS = zeros(dimFS,1);
         coefFS(k) = 1;
         fp.BodyFx = {{'x' 'y'} {fs.fem coefFS}};
         fp.BodyFy = {{'x' 'y'} {[] 0}};
         fem = elModelUpdatePDE(fem,'fp',fp);
         fem = elasticSolvePDE(fem,[]);
         sol{ll+1} = fem.sol;

         fp.BodyFx = {{'x' 'y'} {[] 0}};
         fp.BodyFy = {{'x' 'y'} {fs.fem coefFS}};
         fem = elModelUpdatePDE(fem,'fp',fp);
         fem = elasticSolvePDE(fem,[]);
         sol{ll+2} = fem.sol;
     end
    
    
        nPDE = 2; % two dependent variables, ux and uy
        pb = pde(nPDE);
        % Create a geometry entity
        pg = pdeGeometryFromEdges(dl);
        % Set Diri
    
    
    % 
    
end
    




























